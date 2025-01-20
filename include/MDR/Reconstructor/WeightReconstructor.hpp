#ifndef _MDR_WEIGHT_RECONSTRUCTOR_HPP
#define _MDR_WEIGHT_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "MDR/Decomposer/Decomposer.hpp"
#include "MDR/Interleaver/Interleaver.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/Retriever/Retriever.hpp"
#include "MDR/ErrorEstimator/ErrorEstimator.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/SizeInterpreter/SizeInterpreter.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/RefactorUtils.hpp"
#include <cstdio>
#include "WeightUtils.hpp"
#include "ompSZp_typemanager.h"

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of weight refactor
    template<class T, class Decomposer, class InterleaverT, class InterleaverInt, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class WeightReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        WeightReconstructor(Decomposer decomposer, InterleaverT interleaver, InterleaverInt weight_interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), weight_interleaver(weight_interleaver), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever){}

        T * reconstruct(double tolerance){
            return reconstruct(tolerance, -1);
        }
        // reconstruct data from encoded streams
        T * reconstruct(double tolerance, int max_level=-1){
            // Timer timer;
            // timer.start();
            std::vector<std::vector<double>> level_abs_errors;
            uint8_t target_level = level_error_bounds.size() - 1;
            std::vector<std::vector<double>>& level_errors = level_squared_errors;
            if(std::is_base_of<MaxErrorEstimator<T>, ErrorEstimator>::value){
                // std::cout << "ErrorEstimator is base of MaxErrorEstimator, computing absolute error" << std::endl;
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for(int i=0; i<=target_level; i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }
            else if(std::is_base_of<SquaredErrorEstimator<T>, ErrorEstimator>::value){
                std::cout << "ErrorEstimator is base of SquaredErrorEstimator, using level squared error directly" << std::endl;
            }
            else{
                std::cerr << "Customized error estimator not supported yet" << std::endl;
                exit(-1);
            }
            // timer.end();
            // timer.print("Preprocessing");            

            // timer.start();
            auto prev_level_num_bitplanes(level_num_bitplanes);
            if(max_level == -1 || (max_level >= level_num_bitplanes.size())){
                auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
                // retrieve data
                level_components = retriever.retrieve_level_components(level_sizes, retrieve_sizes, prev_level_num_bitplanes, level_num_bitplanes);                
            }
            else{
                std::vector<std::vector<uint32_t>> tmp_level_sizes;
                std::vector<std::vector<double>> tmp_level_errors;
                std::vector<uint8_t> tmp_level_num_bitplanes;
                for(int i=0; i<=max_level; i++){
                    tmp_level_sizes.push_back(level_sizes[i]);
                    tmp_level_errors.push_back(level_errors[i]);
                    tmp_level_num_bitplanes.push_back(level_num_bitplanes[i]);
                }
                auto retrieve_sizes = interpreter.interpret_retrieve_size(tmp_level_sizes, tmp_level_errors, tolerance, tmp_level_num_bitplanes);
                level_components = retriever.retrieve_level_components(tmp_level_sizes, retrieve_sizes, prev_level_num_bitplanes, tmp_level_num_bitplanes);                
                // add level_num_bitplanes
                for(int i=0; i<=max_level; i++){
                    level_num_bitplanes[i] = tmp_level_num_bitplanes[i];
                }
            }
            // check whether to reconstruct to full resolution
            int skipped_level = 0;
            // for(int i=0; i<=target_level; i++){
            //     if(level_num_bitplanes[target_level - i] != 0){
            //         skipped_level = i;
            //         break;
            //     }
            // }
            // TODO: uncomment skip level to reconstruct low resolution data
            // target_level -= skipped_level;
            // timer.end();
            // timer.print("Interpret and retrieval");
            int reconstruct_level = target_level - skipped_level;
            // std::cout << "skipped_level = " << skipped_level << ", target_level = " << +target_level << std::endl;

            bool success = reconstruct(reconstruct_level, prev_level_num_bitplanes);
            retriever.release();
            if(success){
                current_level = reconstruct_level;
                return data.data();
            }
            else{
                std::cerr << "Reconstruct unsuccessful, return NULL pointer" << std::endl;
                return NULL;
            }
        }

        T * progressive_reconstruct(double tolerance){
            return progressive_reconstruct(tolerance, -1);
        }
        // reconstruct progressively based on available data
        T * progressive_reconstruct(double tolerance, int max_level=-1){
            // std::vector<T> cur_data(data);
            reconstruct(tolerance, max_level);
            // TODO: add resolution changes
            // if(cur_data.size() == data.size()){
            //     for(int i=0; i<data.size(); i++){
            //         data[i] += cur_data[i];
            //     }                
            // }
            // else if(cur_data.size()){
            //     std::cerr << "Reconstruct size changes, not supported yet." << std::endl;
            //     std::cerr << "Sizes before reconstruction: " << cur_data.size() << std::endl;
            //     std::cerr << "Sizes after reconstruction: " << data.size() << std::endl;
            //     exit(0);
            // }
            return data.data();
        }
        // TODO: do not overwrite
        T * recompose_to_full(){
            clear_data(data.data(), current_dimensions, dimensions, dimensions);
            int target_level = level_num.size() - 1;
            std::cout << "recompose to full for " << target_level - current_level << " levels!\n"; 
            std::cout << "dimensions: " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << "\n";
            decomposer.recompose(data.data(), dimensions, target_level - current_level, this->strides); 
            return data.data();
        }

        void load_metadata(){
            uint8_t * metadata = retriever.load_metadata();
            uint8_t const * metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos ++);
            deserialize(metadata_pos, num_dims, dimensions);
            weight_file_size = *reinterpret_cast<const size_t *>(metadata_pos);
            metadata_pos += sizeof(size_t);
            uint8_t num_levels = *(metadata_pos ++);
            deserialize(metadata_pos, num_levels, level_error_bounds);
            // deserialize(metadata_pos, num_levels, level_squared_errors);
            deserialize(metadata_pos, num_levels, level_sizes);
            deserialize(metadata_pos, num_levels, stopping_indices);
            deserialize(metadata_pos, num_levels, level_num);
            negabinary = *(metadata_pos ++);
            level_num_bitplanes = std::vector<uint8_t>(num_levels, 0);
            strides = std::vector<uint32_t>(dimensions.size());
            uint32_t stride = 1;
            for(int i=dimensions.size()-1; i>=0; i--){
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
            free(metadata);
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

        void load_weight(){
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
            memcpy(&bit_count, weight_data_pos, sizeof(unsigned int));
            weight_data_pos += sizeof(unsigned int);
            size_t intArrayLength;
            memcpy(&intArrayLength, weight_data_pos, sizeof(size_t));
            weight_data_pos += sizeof(size_t);
            size_t compressed_weight_size;
            memcpy(&compressed_weight_size, weight_data_pos, sizeof(size_t));
            weight_data_pos += sizeof(size_t);
            compressed_weights.clear();
            compressed_weights.assign(weight_data_pos, weight_data_pos + compressed_weight_size);
            weight_data_pos += compressed_weight_size;
            free(weight_data);
            block_weights.resize(intArrayLength);
        }

        void span_weight(){
            if (compressed_weights.size() != Jiajun_extract_fixed_length_bits(compressed_weights.data(), block_weights.size(), reinterpret_cast<unsigned int*>(block_weights.data()), bit_count)){
                // perror("From WeightedApproximationBasedReconstructor: Error: byteLength != weight_size\n");
            }
            if(dimensions.size() == 1){
                int_weights = fill_block_weight_1D(dimensions[0], block_weights, block_size);
            }
            else if(dimensions.size() == 2){
                int_weights = fill_block_weight_2D(dimensions[0], dimensions[1], block_weights, block_size);
            }
            else{
                int_weights = fill_block_weight_3D(dimensions[0], dimensions[1], dimensions[2], block_weights, block_size);
            }
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
            // for(int i=20153888; i<20153889; i++){
            //     std::cout << "int_weights[" << i << "] = " << int_weights[i] << std::endl;
            // }
            // using level_error_bounds temporarily
            // std::cout << "(int) level_error_bounds.size() = " << (int) level_error_bounds.size() << std::endl;
            propagateWeight(dimensions, (int) level_error_bounds.size() - 1, int_weights);
            write_weight_dat(block_size);
        }

        const std::vector<uint32_t>& get_dimensions(){
            return dimensions;
        }

        const std::vector<uint32_t>& get_current_dimensions(){
            return current_dimensions;
        }

        int get_reconstruct_level(){
            return current_level;
        }

        size_t get_retrieved_size(){
            return retriever.get_retrieved_size();
        }

        size_t get_weight_file_size()
        {
            return weight_file_size;
        }

        int get_max_weight()
        {
            return max_weight;
        }

        std::vector<uint32_t> get_offsets(){
            return retriever.get_offsets();
        }

        std::vector<int> get_int_weights(){
            return int_weights;
        }

        ~WeightReconstructor(){}

        void print() const {
            std::cout << "Weight reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "SizeInterpreter: "; interpreter.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:
        bool reconstruct(uint8_t target_level, const std::vector<uint8_t>& prev_level_num_bitplanes, bool progressive=true){
            auto num_levels = level_num.size();
            auto level_dims = compute_level_dims(dimensions, num_levels - 1);
            auto reconstruct_dimensions = level_dims[target_level];
            
            // print_vec(int_weights);
            // std::cout << "target_level = " << +target_level << ", dims = " << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
            // update with stride
            std::vector<T> cur_data(data);
            memset(data.data(), 0, data.size() * sizeof(T));

            // std::cout << "current_level = " << current_level << std::endl;
            auto level_elements = compute_level_elements(level_dims, target_level);
            std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
            for(int i=0; i<=current_level; i++){
                if(level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                    int * buffer_weight = (int *) malloc(level_elements[i] * sizeof(int));
                    weight_interleaver.interleave(int_weights.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<int*>(buffer_weight));
                    /*for(int k=0; k < level_elements[i]; k++){
                        std::cout << buffer_weight[k] << " ";
                    }
                    std::cout << std::endl;*/
                    int level_exp = 0;
                    if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                    else frexp(level_error_bounds[i], &level_exp);
                    int level_max_weight = compute_max_abs_value(reinterpret_cast<int*>(buffer_weight), level_elements[i]);
                    // std::cout << "level = " << i << ", level_exp = " << level_exp << ", max_weight = " << level_max_weight << std::endl;
                    /*
                    int exp_max_weight = 0;
                    frexp(level_max_weight, &exp_max_weight);
                    // */
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp + level_max_weight, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i, buffer_weight);
                    compressor.decompress_release();
                    interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), this->strides);
                    free(level_decoded_data);                    
                }
            }
            // decompose data to current level
            if(current_level >= 0){
                if(current_level) decomposer.recompose(data.data(), current_dimensions, current_level, this->strides);
                // std::cout << "update data\n";
                // update data with strides
                if(current_dimensions.size() == 1){
                    for(int i=0; i<current_dimensions[0]; i++){
                        data[i] += cur_data[i];
                    }
                }
                else if(current_dimensions.size() == 2){
                    for(int i=0; i<current_dimensions[0]; i++){
                        for(int j=0; j<current_dimensions[1]; j++){
                            data[i*this->strides[0] + j] += cur_data[i*this->strides[0] + j];
                        }
                    }                    
                }
                else if(current_dimensions.size() == 3){
                    for(int i=0; i<current_dimensions[0]; i++){
                        for(int j=0; j<current_dimensions[1]; j++){
                            for(int k=0; k<current_dimensions[2]; k++){
                                data[i*this->strides[0] + j*this->strides[1] + k] += cur_data[i*this->strides[0] + j*this->strides[1] + k];
                            }
                        }
                    }                    
                }
                else{
                    std::cout << "dimension higher than 4 is not supported in update data\n";
                    exit(-1);
                }
            }
            // std::cout << "decompose to target_level\n";
            // decompose data to target level
            for(int i=current_level+1; i<=target_level; i++){
                compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                int * buffer_weight = (int *) malloc(level_elements[i] * sizeof(int));
                weight_interleaver.interleave(int_weights.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<int*>(buffer_weight));
                /*for(int k=0; k < level_elements[i]; k++){
                    std::cout << buffer_weight[k] << " ";
                }
                std::cout << std::endl;*/
                int level_exp = 0;
                if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                else frexp(level_error_bounds[i], &level_exp);
                int level_max_weight = compute_max_abs_value(reinterpret_cast<int*>(buffer_weight), level_elements[i]);
                // std::cout << "level = " << i << ", max_weight = " << level_max_weight << std::endl;
                /*int exp_max_weight = 0;
                frexp(level_max_weight, &exp_max_weight);*/
                auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp + level_max_weight, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i, buffer_weight);
                compressor.decompress_release();
                interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), this->strides);
                free(level_decoded_data);                    
            }
            if(current_level >= 0){
                decomposer.recompose(data.data(), reconstruct_dimensions, target_level - current_level, this->strides);                
            }
            else{
                decomposer.recompose(data.data(), reconstruct_dimensions, target_level, this->strides);
            }
            current_dimensions = reconstruct_dimensions;
            return true;

        }

        void clear_data(T * dst, const std::vector<uint32_t>& coarse_dims, const std::vector<uint32_t>& fine_dims, const std::vector<uint32_t>& dims){
            for(int i=0; i<fine_dims[0]; i++){
                for(int j=0; j<fine_dims[1]; j++){
                    for(int k=0; k<fine_dims[2]; k++){
                        if((i<coarse_dims[0]) && (j<coarse_dims[1]) && (k<coarse_dims[2])){

                        }
                        else{
                            dst[i*dims[1]*dims[2] + j*dims[2] + k] = 0;
                        }
                    }
                }
            }
        }

        Decomposer decomposer;
        InterleaverT interleaver;
        InterleaverInt weight_interleaver;
        Encoder encoder;
        SizeInterpreter interpreter;
        Retriever retriever;
        Compressor compressor;
        size_t weight_file_size = 0;
        unsigned int bit_count;
        std::vector<T> data;
        int block_size;
        int max_weight;
        std::vector<unsigned char> compressed_weights;
        std::vector<int> block_weights;
        std::vector<uint32_t> dimensions;
        std::vector<uint32_t> current_dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> level_num_bitplanes;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<const uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        int current_level = -1;
        std::vector<uint32_t> strides;
        bool negabinary = true;
    public:
        std::vector<int> int_weights;
    };
}
#endif

