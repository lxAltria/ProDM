#ifndef _MDR_FUSE_COMPOSED_RECONSTRUCTOR_NEW_HPP
#define _MDR_FUSE_COMPOSED_RECONSTRUCTOR_NEW_HPP

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

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of composed refactor
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class FuseComposedReconstructor_new : public concepts::ReconstructorInterface<T> {
    public:
        FuseComposedReconstructor_new(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever){}

        T * reconstruct(double tolerance){
            return reconstruct(tolerance, -1);
        }
        // reconstruct data from encoded streams
        T * reconstruct(double tolerance, int max_level=-1){
            // Timer timer;
            // timer.start();
            std::vector<std::vector<double>> level_abs_errors;
            uint8_t target_level = level_error_bounds.size() - 1;

            // std::cout << "level_error_bounds: " << std::endl;
            // for(int i=0; i<level_error_bounds.size(); i++){
            //     std::cout << level_error_bounds[i] << " ";
            // }
            // std::cout << std::endl;
            
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
                // std::cout << "level_num_bitplanes: " << std::endl;
                // for(int i=0; i<level_num_bitplanes.size(); i++){
                //     std::cout << (int)level_num_bitplanes[i] << " ";
                // }
                // std::cout << std::endl;
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
            data = decomposer.recompose(level_buffers, dimensions, target_level - current_level, this->strides); 
            return data.data();
        }

        void load_metadata(){
            uint8_t * metadata = retriever.load_metadata();
            uint8_t const * metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos ++);
            deserialize(metadata_pos, num_dims, dimensions);
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

            MGARD_target_level = (level_error_bounds.size() - 1) / dimensions.size();
            std::cout << "MGARD_target_level = " << MGARD_target_level << std::endl;
            auto level_dims = compute_level_dims_new(dimensions, MGARD_target_level);
            level_buffer_sizes = compute_level_buffers_size(level_dims, MGARD_target_level);

            free(metadata);
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

        std::vector<uint32_t> get_offsets(){
            return retriever.get_offsets();
        }

        ~FuseComposedReconstructor_new(){}

        void print() const {
            std::cout << "New Fuse Composed reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "SizeInterpreter: "; interpreter.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:
        bool reconstruct(uint8_t target_level, const std::vector<uint8_t>& prev_level_num_bitplanes, bool progressive=true){
            auto num_levels = level_num.size();
            auto reconstruct_dimensions = dimensions;
            // update with stride
            std::vector<T> cur_data(data);
            memset(data.data(), 0, data.size() * sizeof(T));
            for(int i=0; i<level_buffer_sizes.size(); i++){
                level_buffers.push_back(std::vector<T>(level_buffer_sizes[i], 0));
            }
            // auto level_elements = compute_level_elements(level_dims, target_level);
            // std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
            for(int i=0; i<=current_level; i++){
                if(level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                    else frexp(level_error_bounds[i], &level_exp);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_buffer_sizes[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    compressor.decompress_release();
                    // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                    // interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), i, target_level, this->strides);
                    // free(level_decoded_data);
                    memcpy(level_buffers[i].data(), level_decoded_data, level_buffer_sizes[i] * sizeof(T));                    
                }
            }
            // decompose data to current level
            if(current_level >= 0){
                if(current_level) data = decomposer.reposition_recompose(level_buffers, current_dimensions, MGARD_target_level, this->strides);
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
            // decompose data to target level
            for(int i=current_level+1; i<=target_level; i++){
                compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                int level_exp = 0;
                if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                else frexp(level_error_bounds[i], &level_exp);
                auto level_decoded_data = encoder.progressive_decode(level_components[i], level_buffer_sizes[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                compressor.decompress_release();
                // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                // interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), i, target_level, this->strides);
                // free(level_decoded_data);          
                memcpy(level_buffers[i].data(), level_decoded_data, level_buffer_sizes[i] * sizeof(T));
            }
            if(current_level >= 0){
                // data = decomposer.reposition_recompose(level_buffers, reconstruct_dimensions, MGARD_target_level, this->strides);      
            }
            else{
                data = decomposer.reposition_recompose(level_buffers, reconstruct_dimensions, MGARD_target_level, this->strides);
            }
            current_dimensions = reconstruct_dimensions;
            level_buffers.clear();
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
        Interleaver interleaver;
        Encoder encoder;
        SizeInterpreter interpreter;
        Retriever retriever;
        Compressor compressor;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<uint32_t> current_dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> level_num_bitplanes;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<const uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<std::vector<T>> level_buffers;
        std::vector<uint32_t> level_buffer_sizes;
        int current_level = -1;
        size_t MGARD_target_level;
        std::vector<uint32_t> strides;
        bool negabinary = true;
    };
}
#endif

