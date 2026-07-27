#ifndef _PDR_COMPOSED_RECONSTRUCTOR_HPP
#define _PDR_COMPOSED_RECONSTRUCTOR_HPP

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

using namespace MDR;

namespace PDR
{
    // an approximation-based scientific data reconstructor: inverse operator of approximation-based refactor
    template <class T, class Approximator, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class ApproximationBasedReconstructor : public concepts::ReconstructorInterface<T>
    {
    public:
        ApproximationBasedReconstructor(Approximator approximator, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
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
            approximator_size = *reinterpret_cast<const size_t *>(metadata_pos);
            metadata_pos += sizeof(size_t);
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

        const std::vector<uint32_t> &get_dimensions()
        {
            return dimensions;
        }

        size_t get_retrieved_size()
        {
            return retriever.get_retrieved_size() + approximator_size;
        }

        std::vector<uint32_t> get_offsets()
        {
            return retriever.get_offsets();
        }

        ~ApproximationBasedReconstructor() {}

        void print() const
        {
            std::cout << "Approximation-based reconstructor with the following components." << std::endl;
            std::cout << "Approximator: ";
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
                std::string approximator_path = retriever.get_directory() + "approximator.dat";
                approximator.reconstruct_approximate(data.data(), dimensions, approximator_path);
                reconstructed = true;
            }
            int i = 0;
            std::vector<uint32_t> level_elements;
            level_elements.push_back(num_elements);
            if (level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0)
            {
                if (mask.empty()){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if (negabinary)
                        frexp(level_error_bounds[i] / 4, &level_exp);
                    else
                        frexp(level_error_bounds[i], &level_exp);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    compressor.decompress_release();
                    for (int i = 0; i < num_elements; i++)
                    {
                        data[i] += level_decoded_data[i];
                    }
                    free(level_decoded_data);
                }
                else{
                    int num_valid_data = std::accumulate(mask.begin(), mask.end(), 0);
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if (negabinary)
                        frexp(level_error_bounds[i] / 4, &level_exp);
                    else
                        frexp(level_error_bounds[i], &level_exp);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], num_valid_data, level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    compressor.decompress_release();
                    int filtered_index = 0;
                    for (int i = 0; i < mask.size(); i++)
                    {
                        if(mask[i]) data[i] += level_decoded_data[filtered_index++];
                        else data[i] = 0;
                    }
                    free(level_decoded_data);
                }
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
        size_t approximator_size = 0;
        std::vector<T> data;
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
        std::vector<unsigned char> mask;
    };
}
#endif
