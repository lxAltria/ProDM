#ifndef _MDR_PARTIAL_ORDERED_COEFFICIENT_PREDICTION_RECONSTRUCTOR_HPP
#define _MDR_PARTIAL_ORDERED_COEFFICIENT_PREDICTION_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "ProDM/Decomposer/MultiLevel/Decomposer.hpp"
#include "ProDM/Decomposer/MultiLevel/Interleaver/Interleaver.hpp"
#include "ProDM/Encoder/BitplaneEncoder.hpp"
#include "ProDM/Retriever/Retriever.hpp"
#include "ProDM/ErrorControl/ErrorEstimator.hpp"
#include "ProDM/ErrorControl/ErrorCollector.hpp"
#include "ProDM/ErrorControl/SizeInterpreter.hpp"
#include "ProDM/Compressor/LevelCompressor.hpp"
#include "ProDM/Utils/RefactorUtils.hpp"

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of composed refactor
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class PartialOrderedCPReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        PartialOrderedCPReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
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
            uint16_t prev_coeff_index = coeff_index;
            if(max_level == -1 || (max_level >= level_num_bitplanes.size())){
                // std::cout << "level_num_bitplanes: " << std::endl;
                // for(int i=0; i<level_num_bitplanes.size(); i++){
                //     std::cout << (int)level_num_bitplanes[i] << " ";
                // }
                // std::cout << std::endl;
                auto retrieve_sizes = interpreter.interpret_retrieve_size_cp_ordered(level_sizes, level_errors, coeff_sizes, combined_coeff_error_perstep, coeff_combined_bps, mdr_target_level, tolerance, level_num_bitplanes, coeff_index);
                // std::cout << "retrieve_sizes: " << std::endl;
                // for(int i=0; i<level_num_bitplanes.size(); i++){
                //     std::cout << retrieve_sizes[i] << " ";
                // }
                // std::cout << std::endl;
                // std::cout << "level retrieve_sizes validation: " << std::endl;
                // for(int i=0; i<level_num_bitplanes.size(); i++){
                //     uint32_t tmp_size = 0;
                //     for(int j=prev_level_num_bitplanes[i]; j<level_num_bitplanes[i]; j++){
                //         tmp_size += level_sizes[i][j];
                //     }
                //     std::cout << tmp_size << " ";
                // }
                // std::cout << std::endl;
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
            // for(size_t i=0; i<level_components.size(); i++){
            //     std::cout << "level_components[" << i << "].size() = " << level_components[i].size() << std::endl;
            // }
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

        void retrieve_metadata(uint8_t* metadata){
            const uint8_t * p = metadata;
            uint8_t num_dims = *(p++);
            deserialize(p, num_dims, dimensions);
            uint8_t num_levels = *(p++);
            deserialize(p, num_levels, level_error_bounds);
            uint8_t interpolation_type = *(p++);
            deserialize(p, num_levels, level_sizes);
            deserialize(p, num_levels, stopping_indices);
            negabinary = (*(p++) != 0);

            uint16_t chunk_num = 0;
            memcpy(&chunk_num, p, sizeof(uint16_t));
            p += sizeof(uint16_t);
            deserialize(p, chunk_num, level_chunk_order);
            deserialize(p, chunk_num, level_error_perstep);

            level_num_bitplanes = std::vector<uint8_t>(num_levels, 0);
            level_num = std::vector<uint32_t>(num_levels, 1);
            level_components = std::vector<std::vector<const uint8_t*>>(num_levels);
            strides = std::vector<uint32_t>(dimensions.size());
            uint32_t stride = 1;
            for(int i=dimensions.size()-1; i>=0; i--){
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
            free(metadata);
        }

        void load_metadata(){
            uint8_t * metadata = retriever.load_metadata();
            const uint8_t * p = metadata;
            uint8_t num_dims = *(p++);
            deserialize(p, num_dims, dimensions);
            uint8_t num_levels = *(p++);
            deserialize(p, num_levels, level_error_bounds);
            uint8_t interpolation_type = *(p++);
            deserialize(p, num_levels, level_sizes);
            deserialize(p, num_levels, level_elements);
            deserialize(p, num_levels, stopping_indices);
            negabinary = (*(p++) != 0);
            deserialize(p, num_dims + 1, decomposed_buffer_dims);
            deserialize(p, num_dims, coeff_interp_directions);
            uint16_t coeff_combined_bp_num = 0;
            memcpy(&coeff_combined_bp_num, p, sizeof(uint16_t));
            p += sizeof(uint16_t);
            deserialize(p, coeff_combined_bp_num, coeff_sizes);
            deserialize(p, coeff_combined_bp_num, combined_coeff_error_perstep);
            deserialize(p, coeff_combined_bp_num, coeff_combined_bps);

            mdr_target_level = num_levels;
            coeff_target_level = 2;
            for(int i=0; i<coeff_interp_directions.size(); i++){
                mdr_target_level -= (coeff_interp_directions[i] == -1) ? 1 : (coeff_target_level + 1);
            }    
            level_retrieve_sizes.resize(num_levels, 0);

            level_num_bitplanes = std::vector<uint8_t>(num_levels, 0);
            level_num = std::vector<uint32_t>(num_levels, 1);
            level_components = std::vector<std::vector<const uint8_t*>>(num_levels);
            strides = std::vector<uint32_t>(dimensions.size());
            uint32_t stride = 1;
            for(int i=dimensions.size()-1; i>=0; i--){
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
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

        size_t get_metadata_size(){
            return retriever.get_metadata_size();
        }

        size_t get_retrieved_size(){
            return retriever.get_retrieved_size();
        }

        std::vector<uint32_t> get_offsets(){
            std::vector<uint32_t> offsets;
            offsets.push_back(retriever.get_offsets());
            return offsets;
        }

        ~PartialOrderedCPReconstructor(){}

        void print() const {
            std::cout << "Ordered coefficient prediction reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "SizeInterpreter: "; interpreter.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:
        bool reconstruct(uint8_t target_level, const std::vector<uint8_t>& prev_level_num_bitplanes, bool progressive=true){
            auto num_levels = level_num.size();
            auto level_dims = compute_level_dims_new(dimensions, 1);
            auto reconstruct_dimensions = level_dims[1];
            // std::cout << "target_level = " << +target_level << ", dims = " << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
            // update with stride
            std::vector<T> cur_data(data);
            memset(data.data(), 0, data.size() * sizeof(T));

            level_buffers.clear();
            level_buffers.resize(num_levels);
            for(int i=0; i<num_levels; i++){
                level_buffers[i].resize(level_elements[i]);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
            }
            decomposed_buffers.clear();
            decomposed_buffers.resize(decomposed_buffer_dims.size());

            for(int i=0; i<=current_level; i++){
                if(level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                    else frexp(level_error_bounds[i], &level_exp);
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    memcpy(level_buffers[i].data(), level_decoded_data, level_elements[i] * sizeof(T));
                    compressor.decompress_release();
                    free(level_decoded_data);                    
                }
                else{
                    memset(level_buffers[i].data(), 0, level_elements[i] * sizeof(T));
                }
            }
            // decompose data to current level
            if(current_level >= 0){
                if(current_level){
                    std::vector<std::vector<T>> tmp_decomposed_buffers;
                    for(int i=0; i<mdr_target_level; i++){
                        tmp_decomposed_buffers.push_back(level_buffers[i]);
                    }
                    auto tmp_data = decomposer.reposition_recompose_split_levels(tmp_decomposed_buffers, decomposed_buffer_dims[0], mdr_target_level - 1);
                    decomposed_buffers[0] = tmp_data;
                    auto coeff_decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
                    uint8_t num_level_1 = mdr_target_level;
                    for(int i=1; i<decomposed_buffer_dims.size(); i++){
                        if(coeff_interp_directions[i - 1] == -1){
                            decomposed_buffers[i] = level_buffers[num_level_1++];
                        }
                        else{
                            coeff_decomposer.direction = coeff_interp_directions[i - 1];
                            std::vector<std::vector<T>> decomposed_coeff_buffers(coeff_target_level + 1);
                            for(int j=0; j<coeff_target_level+1; j++){
                                decomposed_coeff_buffers[j] = level_buffers[num_level_1 + j];
                            }
                            num_level_1 += coeff_target_level + 1;
                            decomposed_buffers[i] = coeff_decomposer.reposition_recompose_split_levels(decomposed_coeff_buffers, decomposed_buffer_dims[i], coeff_target_level);
                        }
                    }
                    data = decomposer.reposition_recompose(decomposed_buffers, current_dimensions, 1, this->strides);
                }
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
                // if(level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                    else frexp(level_error_bounds[i], &level_exp);
                    // std::cout << "level " << i << ", decoding starts" << std::endl;
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    // std::cout << "level " << i << ", decoding done" << std::endl;
                    memcpy(level_buffers[i].data(), level_decoded_data, level_elements[i] * sizeof(T));
                    compressor.decompress_release();
                    free(level_decoded_data);     
                // }
                // else{
                //     memset(level_buffers[i].data(), 0, level_elements[i] * sizeof(T));
                // }        
            }
            if(current_level >= 0){
                // decomposer.recompose(data.data(), reconstruct_dimensions, target_level - current_level, this->strides);                
            }
            else{
                std::vector<std::vector<T>> tmp_decomposed_buffers;
                for(int i=0; i<mdr_target_level; i++){
                    tmp_decomposed_buffers.push_back(level_buffers[i]);
                }
                auto tmp_data = decomposer.reposition_recompose_split_levels(tmp_decomposed_buffers, decomposed_buffer_dims[0], mdr_target_level - 1);
                decomposed_buffers[0] = tmp_data;
                // std::cout << "decomposed_buffers[0].size() = " << decomposed_buffers[0].size() << std::endl;
                auto coeff_decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
                uint8_t num_level_1 = mdr_target_level;
                for(int i=1; i<decomposed_buffer_dims.size(); i++){
                    if(coeff_interp_directions[i - 1] == -1){
                        decomposed_buffers[i] = level_buffers[num_level_1++];
                    }
                    else{
                        coeff_decomposer.direction = coeff_interp_directions[i - 1];
                        // std::cout << "direction = " << (int)coeff_interp_directions[i - 1] << std::endl;
                        std::vector<std::vector<T>> decomposed_coeff_buffers(coeff_target_level + 1);
                        // std::cout << decomposed_buffer_dims[i][0] << " " << decomposed_buffer_dims[i][1] << " " << decomposed_buffer_dims[i][2] << std::endl;
                        for(int j=0; j<coeff_target_level+1; j++){
                            // std::cout << "decomposed_coeff_buffers[" << j << "] = level_buffers[" << num_level_1 + j << "]" << std::endl;
                            decomposed_coeff_buffers[j] = level_buffers[num_level_1 + j];
                        }
                        num_level_1 += coeff_target_level + 1;
                        decomposed_buffers[i] = coeff_decomposer.reposition_recompose_split_levels(decomposed_coeff_buffers, decomposed_buffer_dims[i], coeff_target_level);
                        // std::cout << "decomposed_buffers[" << i << "].size() = " << decomposed_buffers[i].size() << std::endl;
                    }
                }
                // std::cout << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
                data = decomposer.reposition_recompose(decomposed_buffers, reconstruct_dimensions, 1, this->strides);
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
        Interleaver interleaver;
        Encoder encoder;
        SizeInterpreter interpreter;
        Retriever retriever;
        Compressor compressor;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<uint32_t> current_dimensions;
        std::vector<T> level_error_bounds;
        uint16_t num_chunks = 0;
        std::vector<uint8_t> level_num_bitplanes;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<const uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<std::vector<T>> level_buffers;
        std::vector<std::vector<T>> decomposed_buffers;
        std::vector<int> coeff_interp_directions;
        std::vector<std::vector<uint32_t>> decomposed_buffer_dims;
        std::vector<std::vector<uint8_t>> structures;
        std::vector<uint32_t> level_elements;
        int current_level = -1;
        int coeff_target_level;
        int mdr_target_level = 0;
        std::vector<uint32_t> strides;
        bool negabinary = true;
        std::vector<uint8_t> level_chunk_order;
        std::vector<double> level_error_perstep;
        std::vector<uint32_t> chunk_sizes;
        std::vector<size_t> level_retrieve_sizes;
        std::vector<uint32_t> coeff_sizes;
        std::vector<std::vector<uint8_t>> coeff_combined_bps;
        std::vector<double> combined_coeff_error_perstep;
        uint16_t coeff_index = 0;
    };
}
#endif

