#ifndef _PDR_WEIGHTED_COMPOSED_RECONSTRUCTOR_HPP
#define _PDR_WEIGHTED_COMPOSED_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "PDR/Approximator/Approximator.hpp"
#include "MDR/Interleaver/Interleaver.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/Retriever/Retriever.hpp"
#include "MDR/ErrorEstimator/ErrorEstimator.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/SizeInterpreter/SizeInterpreter.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/RefactorUtils.hpp"
#include "WeightUtils.hpp"

using namespace MDR;

namespace PDR
{
    // an approximation-based scientific data reconstructor: inverse operator of approximation-based refactor
    template <class T, class Approximator, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class WeightedApproximationBasedReconstructor : public concepts::ReconstructorInterface<T>
    {
    public:
        WeightedApproximationBasedReconstructor(Approximator approximator, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : approximator(approximator), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever) {}

        T *reconstruct(double tolerance)
        {
            return reconstruct(tolerance, -1);
        }
        // reconstruct data from encoded streams
        T *reconstruct(double tolerance, int max_level = -1)
        {
            // Timer timer;
            // timer.start();
            std::vector<std::vector<double>> level_abs_errors;
            uint8_t target_level = level_error_bounds.size() - 1;
            std::vector<std::vector<double>> &level_errors = level_squared_errors;
            {
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for (int i = 0; i <= target_level; i++)
                {
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }

            // timer.start();
            auto prev_level_num_bitplanes(level_num_bitplanes);
            if (max_level == -1 || (max_level >= level_num_bitplanes.size()))
            {
                auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
                // retrieve data
                level_components = retriever.retrieve_level_components(level_sizes, retrieve_sizes, prev_level_num_bitplanes, level_num_bitplanes);
            }
            else
            {
                std::vector<std::vector<uint32_t>> tmp_level_sizes;
                std::vector<std::vector<double>> tmp_level_errors;
                std::vector<uint8_t> tmp_level_num_bitplanes;
                for (int i = 0; i <= max_level; i++)
                {
                    tmp_level_sizes.push_back(level_sizes[i]);
                    tmp_level_errors.push_back(level_errors[i]);
                    tmp_level_num_bitplanes.push_back(level_num_bitplanes[i]);
                }
                auto retrieve_sizes = interpreter.interpret_retrieve_size(tmp_level_sizes, tmp_level_errors, tolerance, tmp_level_num_bitplanes);
                level_components = retriever.retrieve_level_components(tmp_level_sizes, retrieve_sizes, prev_level_num_bitplanes, tmp_level_num_bitplanes);
                // add level_num_bitplanes
                for (int i = 0; i <= max_level; i++)
                {
                    level_num_bitplanes[i] = tmp_level_num_bitplanes[i];
                }
            }

            bool success = reconstruct(prev_level_num_bitplanes);
            retriever.release();
            if (success)
            {
                return data.data();
            }
            else
            {
                std::cerr << "Reconstruct unsuccessful, return NULL pointer" << std::endl;
                return NULL;
            }
        }

        T *progressive_reconstruct(double tolerance)
        {
            return progressive_reconstruct(tolerance, -1);
        }
        // reconstruct progressively based on available data
        T *progressive_reconstruct(double tolerance, int max_level = -1)
        {
            reconstruct(tolerance, max_level);
            return data.data();
        }

        void load_metadata()
        {
            uint8_t *metadata = retriever.load_metadata();
            uint8_t const *metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos++);
            deserialize(metadata_pos, num_dims, dimensions);
            uint8_t num_levels = *(metadata_pos++);
            deserialize(metadata_pos, num_levels, level_error_bounds);
            deserialize(metadata_pos, num_levels, level_sizes);
            deserialize(metadata_pos, num_levels, stopping_indices);
            deserialize(metadata_pos, num_levels, level_num);
            negabinary = *(metadata_pos++);
            approximator_eb = *reinterpret_cast<const T *>(metadata_pos);
            level_num_bitplanes = std::vector<uint8_t>(num_levels, 0);
            strides = std::vector<uint32_t>(dimensions.size());
            uint32_t stride = 1;
            for (int i = dimensions.size() - 1; i >= 0; i--)
            {
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
            num_elements = stride;
            free(metadata);
        }

        void load_weight()
        {
            std::cout << "Loading Weight" << std::endl;
            string path = retriever.get_directory() + "weight.bin";
            std::cout << "Path: " << path << std::endl;
            FILE *file = fopen(path.c_str(), "r");
            if (file == nullptr)
            {
                perror("Error opening file");
                return;
            }
            fseek(file, 0, SEEK_END);
            uint32_t num_bytes = ftell(file);
            rewind(file);
            uint8_t *weight_data = (uint8_t *)malloc(num_bytes);
            fread(weight_data, 1, num_bytes, file);
            fclose(file);
            uint8_t *weight_data_pos = weight_data;
            memcpy(&block_size, weight_data_pos, sizeof(int));
            weight_data_pos += sizeof(int);
            size_t block_weight_size;
            memcpy(&block_weight_size, weight_data_pos, sizeof(size_t));
            weight_data_pos += sizeof(size_t);
            block_weights.clear();
            block_weights.assign(reinterpret_cast<const int *>(weight_data_pos), reinterpret_cast<const int *>(weight_data_pos) + block_weight_size);
            {
                int max_w = block_weights[0];
                int min_w = block_weights[0];
                for (int i = 1; i < block_weights.size(); i++)
                {
                    if (block_weights[i] > max_w)
                        max_w = block_weights[i];
                    if (block_weights[i] < min_w)
                        min_w = block_weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
            }
            weight_data_pos += block_weight_size * sizeof(int);
            free(weight_data);
        }

        void span_weight()
        {
            if (dimensions.size() == 1)
            {
                int_weights = fill_block_weight_1D(dimensions[0], block_weights, block_size);
            }
            else if (dimensions.size() == 2)
            {
                int_weights = fill_block_weight_2D(dimensions[0], dimensions[1], block_weights, block_size);
            }
            else
            {
                int_weights = fill_block_weight_3D(dimensions[0], dimensions[1], dimensions[2], block_weights, block_size);
            }
            write_weight_dat(block_size);
        }

        void write_weight_dat(const int block_size) const {
            uint32_t weight_size = get_size(int_weights);
            uint8_t * weight_data = (uint8_t *) malloc(weight_size);
            uint8_t * weight_data_pos = weight_data;
            serialize(int_weights, weight_data_pos);
            string path = retriever.get_directory() + "weight_dec.dat";
            std::cout << "Path: " << path << std::endl;
            FILE * file = fopen(path.c_str(), "w");
            if (file == nullptr) {
                perror("Error opening file");
                return;
            }
            fwrite(weight_data, 1, weight_size, file);
            fclose(file);
            free(weight_data);
        }

        const std::vector<uint32_t> &get_dimensions()
        {
            return dimensions;
        }

        size_t get_retrieved_size()
        {
            return retriever.get_retrieved_size();
        }

        std::vector<uint32_t> get_offsets()
        {
            return retriever.get_offsets();
        }

        std::vector<int> get_int_weights()
        {
            return int_weights;
        }

        ~WeightedApproximationBasedReconstructor() {}

        void print() const
        {
            std::cout << "Approximation-based reconstructor with the following components." << std::endl;
            std::cout << "Approximatorr: ";
            approximator.print();
            std::cout << "Encoder: ";
            encoder.print();
            std::cout << "SizeInterpreter: ";
            interpreter.print();
            std::cout << "Retriever: ";
            retriever.print();
        }

    private:
        bool reconstruct(const std::vector<uint8_t> &prev_level_num_bitplanes, bool progressive = true)
        {

            // std::cout << "current_level = " << current_level << std::endl;
            if (!reconstructed)
            {
                // std::string approximator_path = retriever.get_directory() + "approximator.dat";
                // approximator.approximator_file_name = approximator_path;
                approximator.reconstruct_approximate(data.data(), dimensions);
                reconstructed = true;
            }
            int i = 0;
            std::vector<uint32_t> level_elements;
            level_elements.push_back(num_elements);
            if (level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0)
            {
                compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                int level_exp = 0;
                if (negabinary)
                    frexp(level_error_bounds[i] / 4, &level_exp);
                else
                    frexp(level_error_bounds[i], &level_exp);
                int level_max_weight = compute_max_abs_value(int_weights.data(), level_elements[i]);
                auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp + level_max_weight, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i, int_weights.data());
                compressor.decompress_release();
                for (int i = 0; i < num_elements; i++)
                {
                    data[i] += level_decoded_data[i];
                }
                free(level_decoded_data);
            }
            return true;
        }

        Approximator approximator;
        Encoder encoder;
        SizeInterpreter interpreter;
        Retriever retriever;
        Compressor compressor;
        T approximator_eb;
        size_t num_elements;
        std::vector<T> data;
        int block_size;
        std::vector<int> block_weights;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> level_num_bitplanes;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<const uint8_t *>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<uint32_t> strides;
        bool negabinary = true;
        bool reconstructed = false;

    public:
        std::vector<int> int_weights;
    };
}
#endif
