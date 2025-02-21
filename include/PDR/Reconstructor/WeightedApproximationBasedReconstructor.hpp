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
#include "ompSZp_typemanager.c"

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

        T *reconstruct(double tolerance){
            return reconstruct(tolerance, -1);
        }
        // reconstruct data from encoded streams
        T *reconstruct(double tolerance, int max_level = -1){
            // Timer timer;
            // timer.start();
            std::vector<std::vector<double>> level_abs_errors;
            uint8_t target_level = level_error_bounds.size() - 1;
            std::vector<std::vector<double>> &level_errors = level_squared_errors;{
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for (int i = 0; i <= target_level; i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }

            // timer.start();
            auto prev_level_num_bitplanes(level_num_bitplanes);
            if (max_level == -1 || (max_level >= level_num_bitplanes.size())){
                auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
                // retrieve data
                level_components = retriever.retrieve_level_components(level_sizes, retrieve_sizes, prev_level_num_bitplanes, level_num_bitplanes);
            }
            else
            {
                std::vector<std::vector<uint32_t>> tmp_level_sizes;
                std::vector<std::vector<double>> tmp_level_errors;
                std::vector<uint8_t> tmp_level_num_bitplanes;
                for (int i = 0; i <= max_level; i++){
                    tmp_level_sizes.push_back(level_sizes[i]);
                    tmp_level_errors.push_back(level_errors[i]);
                    tmp_level_num_bitplanes.push_back(level_num_bitplanes[i]);
                }
                auto retrieve_sizes = interpreter.interpret_retrieve_size(tmp_level_sizes, tmp_level_errors, tolerance, tmp_level_num_bitplanes);
                level_components = retriever.retrieve_level_components(tmp_level_sizes, retrieve_sizes, prev_level_num_bitplanes, tmp_level_num_bitplanes);
                // add level_num_bitplanes
                for (int i = 0; i <= max_level; i++){
                    level_num_bitplanes[i] = tmp_level_num_bitplanes[i];
                }
            }

            bool success = reconstruct(prev_level_num_bitplanes);
            retriever.release();
            if (success){
                return data.data();
            }
            else{
                std::cerr << "Reconstruct unsuccessful, return NULL pointer" << std::endl;
                return NULL;
            }
        }

        T *progressive_reconstruct(double tolerance){
            return progressive_reconstruct(tolerance, -1);
        }
        // reconstruct progressively based on available data
        T *progressive_reconstruct(double tolerance, int max_level = -1){
            reconstruct(tolerance, max_level);
            return data.data();
        }

        void load_metadata(){
            uint8_t *metadata = retriever.load_metadata();
            uint8_t const *metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos++);
            deserialize(metadata_pos, num_dims, dimensions);
            approximator_size = *reinterpret_cast<const size_t *>(metadata_pos);
            metadata_pos += sizeof(size_t);
            weight_file_size = *reinterpret_cast<const size_t *>(metadata_pos);
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
            for (int i = dimensions.size() - 1; i >= 0; i--){
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
            num_elements = stride;
            free(metadata);
        }

        void load_weight(){
            std::cout << "Loading Weight" << std::endl;
            string path = retriever.get_directory() + "weight.bin";
            std::cout << "Path: " << path << std::endl;
            FILE *file = fopen(path.c_str(), "r");
            if (file == nullptr){
                perror("Error opening file\n");
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
            memcpy(&bit_count, weight_data_pos, sizeof(unsigned int));
            weight_data_pos += sizeof(unsigned int);
            size_t intArrayLength;
            memcpy(&intArrayLength, weight_data_pos, sizeof(size_t));
            weight_data_pos += sizeof(size_t);
            memcpy(&ZSTD_weight_size, weight_data_pos, sizeof(uint32_t));
            weight_data_pos += sizeof(uint32_t);
            ZSTD_weights = (uint8_t*)malloc(ZSTD_weight_size);
            memcpy(ZSTD_weights, weight_data_pos, ZSTD_weight_size);
            weight_data_pos += ZSTD_weight_size;
            free(weight_data);
            block_weights.resize(intArrayLength);
        }

        void span_weight(){
            size_t byteLength = 0;
            unsigned int byte_count = bit_count / 8;
            unsigned int remainder_bit = bit_count % 8;
            if (remainder_bit == 0){
                byteLength = byte_count * block_weights.size() + 1;
            }
            else{
                size_t tmp = remainder_bit * block_weights.size();
                byteLength = byte_count * block_weights.size() + (tmp - 1) / 8 + 1;
            }
            compressed_weights.clear();
            compressed_weights.resize(byteLength);
            uint8_t* tmp_data = nullptr;
            uint32_t tmp_size = ZSTD::decompress(ZSTD_weights, ZSTD_weight_size, &tmp_data);

            if (tmp_data != nullptr) {
                compressed_weights.assign(tmp_data, tmp_data + tmp_size);
                free(tmp_data);
            }

            if (compressed_weights.size() != Jiajun_extract_fixed_length_bits(compressed_weights.data(), block_weights.size(), reinterpret_cast<unsigned int*>(block_weights.data()), bit_count)){
                // perror("From WeightedApproximationBasedReconstructor: Error: byteLength != weight_size\n");
            }
            if (dimensions.size() == 1){
                int_weights = fill_block_weight_1D(dimensions[0], block_weights, block_size);
            }
            else if (dimensions.size() == 2){
                int_weights = fill_block_weight_2D(dimensions[0], dimensions[1], block_weights, block_size);
            }
            else{
                int_weights = fill_block_weight_3D(dimensions[0], dimensions[1], dimensions[2], block_weights, block_size);
            }
            // write_weight_dat(block_size);
            {
                std::cout << "Spanning Weights" << std::endl;
                int max_w = int_weights[0];
                int min_w = int_weights[0];
                for(int i=1; i< int_weights.size(); i++){
                    if(int_weights[i] > max_w) max_w = int_weights[i];
                    if(int_weights[i] < min_w) min_w = int_weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
                max_weight = max_w;
            }
            // write_weight_dat(block_size);
            free(ZSTD_weights);
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

        const std::vector<uint32_t> &get_dimensions(){
            return dimensions;
        }

        size_t get_retrieved_size(){
            return retriever.get_retrieved_size() + approximator_size;
        }

        size_t get_weight_file_size(){
            return weight_file_size;
        }

        int get_max_weight(){
            return max_weight;
        }

        std::vector<uint32_t> get_offsets(){
            return retriever.get_offsets();
        }

        std::vector<int> get_int_weights(){
            return int_weights;
        }

        ~WeightedApproximationBasedReconstructor() {}

        void print() const{
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
            if (!reconstructed){
                std::string approximator_path = retriever.get_directory() + "approximator.dat";
                approximator.reconstruct_approximate(data.data(), dimensions, approximator_path);
                reconstructed = true;
            }
            int i = 0;
            std::vector<uint32_t> level_elements;
            level_elements.push_back(num_elements);
            if (level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                if (mask.empty()){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if (negabinary)
                        frexp(level_error_bounds[i] / 4, &level_exp);
                    else
                        frexp(level_error_bounds[i], &level_exp);
                    int level_max_weight = compute_max_abs_value(int_weights.data(), level_elements[i]);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp + level_max_weight, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i, int_weights.data());
                    compressor.decompress_release();
                    for (int i = 0; i < num_elements; i++){
                        data[i] += level_decoded_data[i];
                    }
                    free(level_decoded_data);
                }
                else{
                    int num_valid_data = std::accumulate(mask.begin(), mask.end(), 0);
                    std::vector<int> filtered_int_weights(num_valid_data);
                    int filtered_index = 0;
                    for(int i=0; i<mask.size(); i++){
                        if(mask[i]){
                            filtered_int_weights[filtered_index++]=int_weights[i];
                        }
                    }
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if (negabinary)
                        frexp(level_error_bounds[i] / 4, &level_exp);
                    else
                        frexp(level_error_bounds[i], &level_exp);
                    int level_max_weight = compute_max_abs_value(filtered_int_weights.data(), num_valid_data);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], num_valid_data, level_exp + level_max_weight, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i, filtered_int_weights.data());
                    compressor.decompress_release();
                    filtered_index = 0;
                    for (int i = 0; i < mask.size(); i++){
                        if (mask[i]) data[i] += level_decoded_data[filtered_index++];
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
        size_t weight_file_size = 0;
        uint8_t* ZSTD_weights;
        uint32_t ZSTD_weight_size;
        unsigned int bit_count;
        std::vector<T> data;
        int block_size;
        int max_weight;
        std::vector<unsigned char> compressed_weights;
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
        std::vector<unsigned char> mask;
    };
}
#endif
