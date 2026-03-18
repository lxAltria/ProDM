#ifndef _MDR_ORDERED_RECONSTRUCTOR_NEW_HPP
#define _MDR_ORDERED_RECONSTRUCTOR_NEW_HPP

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
    class OrderedReconstructor_new : public concepts::ReconstructorInterface<T> {
    public:
        OrderedReconstructor_new(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever){}

        T * reconstruct_from_buffer(double tolerance, const uint8_t *buffer)
        {
            if (!buffer_initialized) {
                if (buffer == NULL) {
                    std::cerr << "reconstruct_from_buffer: buffer is NULL" << std::endl;
                    return NULL;
                }
                buffer_base = buffer;
                const uint8_t *p = buffer_base;

                // get metadata_size
                std::memcpy(&metadata_size_from_buffer, p, sizeof(uint32_t));
                p += sizeof(uint32_t);

                const uint8_t *metadata = p;
                data_base_from_buffer = p + metadata_size_from_buffer;

                // ---- get metadata and initialize the structure ----
                {
                    const uint8_t *mp = metadata;

                    uint8_t num_dims = *(mp++);
                    deserialize(mp, num_dims, dimensions);

                    uint8_t num_levels = *(mp++);
                    deserialize(mp, num_levels, level_error_bounds);
                    deserialize(mp, num_levels, level_sizes);
                    deserialize(mp, num_levels, stopping_indices);

                    negabinary = (*(mp++) != 0);

                    uint16_t chunk_num = 0;
                    std::memcpy(&chunk_num, mp, sizeof(uint16_t));
                    mp += sizeof(uint16_t);
                    deserialize(mp, chunk_num, chunk_order);
                    deserialize(mp, chunk_num, error_perstep);

                    // progressive
                    level_num_bitplanes = std::vector<uint8_t>(num_levels, 0);
                    level_num           = std::vector<uint32_t>(num_levels, 1);
                    level_components    =
                        std::vector<std::vector<const uint8_t *>>(num_levels);

                    strides = std::vector<uint32_t>(dimensions.size());
                    uint32_t stride = 1;
                    for (int i = static_cast<int>(dimensions.size()) - 1; i >= 0; i--) {
                        strides[i] = stride;
                        stride *= dimensions[i];
                    }
                    data = std::vector<T>(stride, 0);

                    num_chunks = 0;
                    chunk_sizes.clear();
                    current_level = -1;
                }

                buffer_initialized = true;
            } else {
                if (buffer != NULL && buffer != buffer_base) {
                    std::cerr << "reconstruct_from_buffer: buffer pointer changed, "
                                "progressive reconstruction assumes same buffer!"
                            << std::endl;
                }
            }

            if (!error_preprocessed) {
                std::vector<std::vector<double>> level_abs_errors;
                uint8_t target_level =
                    static_cast<uint8_t>(level_error_bounds.size() - 1);
                std::vector<std::vector<double>> &level_errors =
                    level_squared_errors;

                if (std::is_base_of<MaxErrorEstimator<T>, ErrorEstimator>::value) {
                    std::cout << "Using absolute error" << std::endl;
                    MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                    level_abs_errors.clear();
                    for (int i = 0; i <= target_level; i++) {
                        auto collected_error =
                            collector.collect_level_error(NULL, 0,
                                                        level_sizes[i].size(),
                                                        level_error_bounds[i]);
                        level_abs_errors.push_back(collected_error);
                    }
                    level_errors = level_abs_errors;
                }
                else if (std::is_base_of<SquaredErrorEstimator<T>, ErrorEstimator>::value) {
                    std::cout << "Using squared error" << std::endl;
                }
                else {
                    std::cerr << "Customized error estimator not supported yet"
                            << std::endl;
                    exit(-1);
                }

                error_preprocessed = true;
            }

            // ---------- add chunks based on tolerance ----------
            auto prev_level_num_bitplanes(level_num_bitplanes);
            size_t retrieve_size = 0;
            size_t prev_num_chunks = num_chunks; // 0 for non-progressive
            double best_error = error_perstep.back();
            if (tolerance < best_error) {
                tolerance = best_error;
            }
            if (prev_num_chunks > 0 && error_perstep[prev_num_chunks - 1] <= tolerance) {
                return data.data();
            }

            for (size_t i = prev_num_chunks; i < chunk_order.size(); i++)
            {
                size_t lv =
                    chunk_order[i];
                size_t sz =
                    level_sizes[lv][level_num_bitplanes[lv]++];
                chunk_sizes.push_back(static_cast<uint32_t>(sz));
                retrieve_size += sz;
                if (error_perstep[i] <= tolerance)
                {
                    num_chunks = static_cast<uint16_t>(i + 1);
                    break;
                }
            }
            
            if (retrieve_size > 0 && num_chunks == prev_num_chunks) {
                num_chunks = static_cast<uint16_t>(chunk_order.size());
            }

            const uint8_t *ordered_components = data_base_from_buffer;

            // skip used chunks 
            size_t offset = 0;
            for (size_t i = 0; i < prev_num_chunks; ++i) {
                offset += chunk_sizes[i];
            }

            // only need to save newly added chunks
            level_components.clear();
            level_components =
                std::vector<std::vector<const uint8_t *>>(level_num.size());

            for (size_t i = prev_num_chunks; i < num_chunks; i++) {
                size_t lv = chunk_order[i];
                level_components[lv].push_back(ordered_components + offset);
                offset += chunk_sizes[i];
            }

            uint8_t target_level =
                static_cast<uint8_t>(level_error_bounds.size() - 1);
            int skipped_level = 0;
            for (int i = 0; i <= target_level; i++) {
                if (level_num_bitplanes[target_level - i] != 0) {
                    skipped_level = i;
                    break;
                }
            }
            int reconstruct_level = static_cast<int>(target_level) - skipped_level;

            bool success =
                reconstruct(static_cast<uint8_t>(reconstruct_level),
                            prev_level_num_bitplanes);

            if (success) {
                current_level = reconstruct_level;
                return data.data();
            } else {
                std::cerr << "Reconstruct unsuccessful, return NULL pointer"
                        << std::endl;
                return NULL;
            }
        }


        // reconstruct data from encoded streams
        T * reconstruct(double tolerance){
            // Timer timer;
            // timer.start();
            std::vector<std::vector<double>> level_abs_errors;
            // uint8_t target_level = level_error_bounds.size() - 1;
            uint8_t target_level = (level_error_bounds.size() - 1) / dimensions.size();
            std::vector<std::vector<double>>& level_errors = level_squared_errors;
            if(std::is_base_of<MaxErrorEstimator<T>, ErrorEstimator>::value){
                std::cout << "Using absolute error" << std::endl;
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for(int i=0; i<level_sizes.size(); i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }
            else if(std::is_base_of<SquaredErrorEstimator<T>, ErrorEstimator>::value){
                std::cout << "Using squared error" << std::endl;
            }
            else{
                std::cerr << "Customized error estimator not supported yet" << std::endl;
                exit(-1);
            }
            // timer.end();
            // timer.print("Preprocessing");            
            // timer.start();

            auto prev_level_num_bitplanes(level_num_bitplanes);
            size_t retrieve_size = 0;
            size_t prev_num_chunks = num_chunks;            
            double best_error = error_perstep.back();
            if (tolerance < best_error) {
                tolerance = best_error;
            }
            if (prev_num_chunks > 0 && error_perstep[prev_num_chunks - 1] <= tolerance) {
                return data.data();
            }

            double estimated_error = 0;
            for (size_t i = prev_num_chunks; i < chunk_order.size(); i++)
            {
                size_t lv = chunk_order[i];
                size_t sz = level_sizes[lv][level_num_bitplanes[lv]++];
                chunk_sizes.push_back(sz);
                retrieve_size += sz;
                if (error_perstep[i] <= tolerance) {
                    estimated_error = error_perstep[i];
                    num_chunks = i + 1;
                    break;
                }
            }
            
            std::cout << "****Tolerance = " << tolerance << ", estimated_error = " << estimated_error << ", num_chunks = " << num_chunks - 1 << ", level_num_bitplanes" << std::endl;
            for(size_t i=0; i<level_num_bitplanes.size(); i++){
                std::cout << (int)level_num_bitplanes[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "level_errors " << std::endl;
            for(size_t i=0; i<level_num_bitplanes.size(); i++){
                std::cout << level_errors[i][level_num_bitplanes[i]] << " ";
            }
            std::cout << std::endl;
            if (retrieve_size > 0 && num_chunks == prev_num_chunks) {
                num_chunks = chunk_order.size();
            }
            uint8_t* ordered_components = retriever.retrieve_components(retrieve_size);
            size_t offset = 0;

            level_components.clear();
            level_components = std::vector<std::vector<const uint8_t*>>(level_num.size());
            
            for (size_t i = prev_num_chunks; i < num_chunks; i++)
            {
                size_t lv = chunk_order[i];
                level_components[lv].push_back(ordered_components + offset);
                offset += chunk_sizes[i];
            }
                        
            // check whether to reconstruct to full resolution
            int skipped_level = 0;
            for(int i=0; i<=target_level; i++){
                if(level_num_bitplanes[target_level - i] != 0){
                    skipped_level = i;
                    break;
                }
            }
            // TODO: uncomment skip level to reconstruct low resolution data
            // target_level -= skipped_level;
            // timer.end();
            // timer.print("Interpret and retrieval");
            int reconstruct_level = level_error_bounds.size() - 1;
            // std::cout << "skipped_level = " << skipped_level << ", target_level = " << +target_level << std::endl;

            bool success = reconstruct(target_level, prev_level_num_bitplanes);
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
            reconstruct(tolerance);
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
            deserialize(p, num_levels, level_sizes);
            deserialize(p, num_levels, stopping_indices);
            negabinary = (*(p++) != 0);

            uint16_t chunk_num = 0;
            memcpy(&chunk_num, p, sizeof(uint16_t));
            p += sizeof(uint16_t);
            deserialize(p, chunk_num, chunk_order);
            deserialize(p, chunk_num, error_perstep);

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
            deserialize(p, num_levels, level_sizes);
            deserialize(p, num_levels, level_elements);
            deserialize(p, num_levels, stopping_indices);
            negabinary = (*(p++) != 0);

            uint16_t chunk_num = 0;
            memcpy(&chunk_num, p, sizeof(uint16_t));
            p += sizeof(uint16_t);
            deserialize(p, chunk_num, chunk_order);
            deserialize(p, chunk_num, error_perstep);

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

        size_t get_retrieved_size(){
            return retriever.get_retrieved_size();
        }

        std::vector<uint32_t> get_offsets(){
            return retriever.get_offsets();
        }

        ~OrderedReconstructor_new(){}

        void print() const {
            std::cout << "Composed reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "SizeInterpreter: "; interpreter.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:        
        bool reconstruct(uint8_t target_level, const std::vector<uint8_t>& prev_level_num_bitplanes, bool progressive=true){
            auto num_levels = level_num.size();
            auto level_dims = compute_level_dims(dimensions, target_level);
            auto reconstruct_dimensions = level_dims[target_level];
            // std::cout << "target_level = " << +target_level << ", dims = " << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
            // update with stride
            std::vector<T> cur_data(data);
            memset(data.data(), 0, data.size() * sizeof(T));
            // std::cout << "Test 1" << std::endl;
            level_buffers.resize(num_levels);
            for(int i=0; i<num_levels; i++){
                level_buffers[i].resize(level_elements[i]);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
            }

            // std::cout << "current_level = " << current_level << std::endl;
            // auto level_elements = compute_level_elements(level_dims, target_level);
            // std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
            // std::cout << "Test 2" << std::endl;
            // std::cout << "current level =" << (int)current_level << std::endl;
            for(int i=0; i<=current_level; i++){
                // std::cout << "i=" << i << " ";
                if(level_num_bitplanes[i] - prev_level_num_bitplanes[i] > 0){
                    // std::cout << "Test 2b" << " ";
                    compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                    int level_exp = 0;
                    if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                    else frexp(level_error_bounds[i], &level_exp);
                    // std::cout << "Test 2c" << " ";
                    auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                    compressor.decompress_release();
                    memcpy(level_buffers[i].data(), level_decoded_data, level_elements[i] * sizeof(T));
                    // std::cout << "Test 2d" << std::endl;
                    // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                    // interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), this->strides);
                    free(level_decoded_data);                    
                }
                else{
                    memset(level_buffers[i].data(), 0, level_elements[i] * sizeof(T));
                }
            }
            // std::cout << "Test 3" << std::endl;
            // decompose data to current level
            if(current_level >= 0){
                if(current_level) {
                    data = decomposer.reposition_recompose(level_buffers, reconstruct_dimensions, target_level, this->strides);
                }
                // std::cout << "update data\n";
                // update data with strides
                if(current_dimensions.size() == 1){
                    for(int i=0; i<current_dimensions[0]; i++){
                        data[i] += cur_data[i];
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
            // std::cout << "Test 4" << std::endl;
            // std::cout << "target level =" << (int)target_level << std::endl;
            // std::cout << "decompose to target_level\n";
            // decompose data to target level
            for(int i=current_level+1; i<num_levels; i++){
                // std::cout << "i=" << i << " ";
                compressor.decompress_level(level_components[i], level_sizes[i], prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], stopping_indices[i]);
                int level_exp = 0;
                if(negabinary) frexp(level_error_bounds[i] / 4, &level_exp);
                else frexp(level_error_bounds[i], &level_exp);
                auto level_decoded_data = encoder.progressive_decode(level_components[i], level_elements[i], level_exp, prev_level_num_bitplanes[i], level_num_bitplanes[i] - prev_level_num_bitplanes[i], i);
                compressor.decompress_release();
                memcpy(level_buffers[i].data(), level_decoded_data, level_elements[i] * sizeof(T));
                // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                // interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data(), this->strides);
                free(level_decoded_data);                    
            }
            // std::cout << "Test 5" << std::endl;
            if(current_level >= 0){
                // decomposer.recompose(data.data(), reconstruct_dimensions, target_level - current_level, this->strides);                
            }
            else{
                data = decomposer.reposition_recompose(level_buffers, reconstruct_dimensions, target_level, this->strides);
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
        int current_level = -1;
        std::vector<uint32_t> strides;
        bool negabinary = true;
        std::vector<uint8_t> chunk_order;
        std::vector<double> error_perstep;
        std::vector<uint32_t> chunk_sizes;
        std::vector<uint32_t> level_elements;
        std::vector<std::vector<T>> level_buffers;

        bool buffer_initialized = false;
        bool error_preprocessed = false;
        const uint8_t *buffer_base = nullptr;
        const uint8_t *data_base_from_buffer = nullptr;
        uint32_t metadata_size_from_buffer = 0;
    };
}
#endif