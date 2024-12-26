#ifndef _MDR_WEIGHT_REFACTOR_HPP
#define _MDR_WEIGHT_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "MDR/Decomposer/Decomposer.hpp"
#include "MDR/Interleaver/Interleaver.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/Writer/Writer.hpp"
#include "MDR/RefactorUtils.hpp"
#include "WeightUtils.hpp"
#include <cstdio>

namespace MDR {
    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, and error collector
    template<class T, class Decomposer, class InterleaverT, class InterleaverInt, class Encoder, class Compressor, class ErrorCollector, class Writer>
    class WeightRefactor : public concepts::RefactorInterface<T> {
    public:
        WeightRefactor(Decomposer decomposer, InterleaverT interleaver, InterleaverInt weight_interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), weight_interleaver(weight_interleaver), encoder(encoder), compressor(compressor), collector(collector), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes, int max_weight = 4, const int block_size = 1){
            Timer timer;
            timer.start();
            dimensions = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            weights = data;
            // if refactor successfully
            if(refactor(target_level, num_bitplanes, max_weight, block_size)){
                timer.end();
                timer.print("Refactor");
                timer.start();
                level_num = writer.write_level_components(level_components, level_sizes);
                timer.end();
                timer.print("Write");                
            }

            write_metadata();
            write_weight(block_size);
            for(int i=0; i<level_components.size(); i++){
                for(int j=0; j<level_components[i].size(); j++){
                    free(level_components[i][j]);                    
                }
            }
        }

        void write_metadata() const {
            uint32_t metadata_size = sizeof(uint8_t) + get_size(dimensions) // dimensions
                            + sizeof(uint8_t) + get_size(level_error_bounds) 
                            // + get_size(level_squared_errors) 
                            + get_size(level_sizes) // level information
                            + get_size(stopping_indices) + get_size(level_num) + 1; // one byte for whether negabinary encoding is used 
            uint8_t * metadata = (uint8_t *) malloc(metadata_size);
            uint8_t * metadata_pos = metadata;
            *(metadata_pos ++) = (uint8_t) dimensions.size();
            serialize(dimensions, metadata_pos);
            *(metadata_pos ++) = (uint8_t) level_error_bounds.size();
            serialize(level_error_bounds, metadata_pos);
            // serialize(level_squared_errors, metadata_pos);
            serialize(level_sizes, metadata_pos);
            serialize(stopping_indices, metadata_pos);
            serialize(level_num, metadata_pos);
            *(metadata_pos ++) = (uint8_t) negabinary;
            writer.write_metadata(metadata, metadata_size);
            free(metadata);
        }

        void write_weight(const int block_size) const {
            {
                int max_w = block_weights[0];
                int min_w = block_weights[0];
                for(int i=1; i< block_weights.size(); i++){
                    if(block_weights[i] > max_w) max_w = block_weights[i];
                    if(block_weights[i] < min_w) min_w = block_weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
            }
 
            uint32_t weight_size = sizeof(int) + sizeof(size_t) + get_size(block_weights);
            uint8_t * weight_data = (uint8_t *) malloc(weight_size);
            uint8_t * weight_data_pos = weight_data;
            memcpy(weight_data_pos, &block_size, sizeof(int));
            weight_data_pos += sizeof(int);
            size_t block_weight_size = block_weights.size();
            memcpy(weight_data_pos, &block_weight_size, sizeof(size_t));  
            weight_data_pos += sizeof(size_t);
            serialize(block_weights, weight_data_pos);
            string path = writer.get_directory() + "/weight.bin";
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

        void write_weight_dat(const int block_size) const {
            uint32_t weight_size = get_size(int_weights);
            uint8_t * weight_data = (uint8_t *) malloc(weight_size);
            uint8_t * weight_data_pos = weight_data;
            serialize(int_weights, weight_data_pos);
            string path = writer.get_directory() + "weight.dat";
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

        ~WeightRefactor(){}

        void print() const {
            std::cout << "Weight refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Weight_interleaver: "; weight_interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        bool refactor(uint8_t target_level, uint8_t num_bitplanes, int max_weight, const int block_size){
            uint8_t max_level = log2(*min_element(dimensions.begin(), dimensions.end())) - 1;
            if(target_level > max_level){
                std::cerr << "Target level is higher than " << max_level << std::endl;
                return false;
            }
            // Timer timer;
            // decompose data hierarchically
            // timer.start();
            if (dimensions.size() == 1){
                std::cout << "assign_block_value_1D\n";
                assign_block_value_1D(dimensions[0], block_size, weights.data());
            }
            else if (dimensions.size() == 2){
                assign_block_value_2D(dimensions[0], dimensions[1], dimensions[1], block_size, weights.data());
            }
            else if (dimensions.size() == 3){
                assign_block_value_3D(dimensions[0], dimensions[1], dimensions[2], dimensions[1]*dimensions[2], dimensions[2], block_size, weights.data());
            }
            {
                int max_w = weights[0];
                int min_w = weights[0];
                for(int i=1; i< weights.size(); i++){
                    if(weights[i] > max_w) max_w = weights[i];
                    if(weights[i] < min_w) min_w = weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
            }
            int_weights = normalize_weights(weights, max_weight);
            {
                std::cout << "Normalizing Weights" << std::endl;
                int max_w = int_weights[0];
                int min_w = int_weights[0];
                for(int i=1; i< int_weights.size(); i++){
                    if(int_weights[i] > max_w) max_w = int_weights[i];
                    if(int_weights[i] < min_w) min_w = int_weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
            }
            if(dimensions.size() == 1){
                block_weights = get_block_weight_1D(dimensions[0], int_weights, block_size);
            }
            else if(dimensions.size() == 2){
                block_weights = get_block_weight_2D(dimensions[0], dimensions[1], int_weights, block_size);
            }
            write_weight_dat(block_size);
            propagateWeight(dimensions, (int) target_level, int_weights);
            {
                int max_w = int_weights[0];
                int min_w = int_weights[0];
                for(int i=1; i< weights.size(); i++){
                    if(int_weights[i] > max_w) max_w = int_weights[i];
                    if(int_weights[i] < min_w) min_w = int_weights[i];
                }
                std::cout << min_w << " " << max_w << std::endl;
            }
            // exit(-1);
            // print_vec(int_weights);
            decomposer.decompose(data.data(), dimensions, target_level);
            // MGARD::writefile("decomposed_coeff.dat", data.data(), data.size());
            // timer.end();
            // timer.print("Decompose");

            // encode level by level
            level_error_bounds.clear();
            level_squared_errors.clear();
            level_components.clear();
            level_sizes.clear();
            auto level_dims = compute_level_dims(dimensions, target_level);
            auto level_elements = compute_level_elements(level_dims, target_level);
            std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
            SquaredErrorCollector<T> s_collector = SquaredErrorCollector<T>();
            size_t compressed_size = 0;
            for(int i=0; i<=target_level; i++){
                // timer.start();
                const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
                int * buffer_weight = (int *) malloc(level_elements[i] * sizeof(int));
                // extract level i component
                interleaver.interleave(data.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
                weight_interleaver.interleave(int_weights.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<int*>(buffer_weight));
                /*for(int k=0; k < level_elements[i]; k++){
                    std::cout << buffer_weight[k] << " ";
                }
                std::cout << std::endl;*/
                // compute max coefficient as level error bound
                T level_max_error = compute_max_abs_value(reinterpret_cast<T*>(buffer), level_elements[i]);
                // /*
                int level_max_weight = compute_max_abs_value(reinterpret_cast<int*>(buffer_weight), level_elements[i]);
                // */
                // std::cout << "\nlevel " << i << " max error = " << level_max_error << std::endl;
                // MGARD::writefile(("level_" + std::to_string(i) + "_coeff.dat").c_str(), buffer, level_elements[i]);
                if(negabinary) level_error_bounds.push_back(level_max_error * 4);
                else level_error_bounds.push_back(level_max_error);
                // timer.end();
                // timer.print("Interleave");
                // collect errors
                // auto collected_error = s_collector.collect_level_error(buffer, level_elements[i], num_bitplanes, level_max_error);
                // level_squared_errors.push_back(collected_error);
                // encode level data
                // timer.start();
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                /*
                int exp_max_weight = 0;
                frexp(level_max_weight, &exp_max_weight);
                // */
                std::vector<uint32_t> stream_sizes;
                //std::vector<double> level_sq_err;
                auto streams = encoder.encode(buffer, level_elements[i], level_exp + level_max_weight, num_bitplanes, stream_sizes, buffer_weight);
                free(buffer);
                // level_squared_errors.push_back(level_sq_err);
                // timer.end();
                // timer.print("Encoding");
                // timer.start();
                // lossless compression
                uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                stopping_indices.push_back(stopping_index);
                // record encoded level data and size
                level_components.push_back(streams);
                level_sizes.push_back(stream_sizes);
                // timer.end();
                // timer.print("Lossless time");
            }
            // print_vec("level sizes", level_sizes);
            return true;
        }

        Decomposer decomposer;
        InterleaverT interleaver;
        InterleaverInt weight_interleaver;
        Encoder encoder;
        Compressor compressor;
        ErrorCollector collector;
        Writer writer;
        std::vector<T> data;
        std::vector<T> weights;
        std::vector<int> int_weights;
        std::vector<int> block_weights;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
    public:
        bool negabinary = false;
    };
}
#endif

