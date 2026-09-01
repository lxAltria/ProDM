#ifndef _MDR_COEFFICIENT_PREDICTION_REFACTOR_HPP
#define _MDR_COEFFICIENT_PREDICTION_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "ProDM/Decomposer/MultiLevel/Decomposer.hpp"
#include "ProDM/Decomposer/MultiLevel/Interleaver/Interleaver.hpp"
#include "ProDM/Encoder/BitplaneEncoder.hpp"
#include "ProDM/ErrorControl/Collector/ErrorCollector.hpp"
#include "ProDM/Compressor/LevelCompressor.hpp"
#include "ProDM/Writer/Writer.hpp"
#include "ProDM/Utils/RefactorUtils.hpp"
#include "ProDM/Decomposer/MultiLevel/Tuner/Tuner.hpp"

namespace MDR {
    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, and error collector
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
    class CPRefactor : public concepts::RefactorInterface<T> {
    public:
        CPRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), collector(collector), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){
            Timer timer;
            timer.start();
            dimensions = dims;
            size_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            if(num_elements > (size_t) UINT32_MAX){
                std::cerr << "Data with more than 2^32 elements is not supported." << std::endl;
                exit(-1);
            }
            data = std::vector<T>(data_, data_ + num_elements);
            // if refactor successfully
            if(refactor(target_level, num_bitplanes)){
                // timer.end();
                // timer.print("Refactor");
                // timer.start();
                level_num = writer.write_level_components(level_components, level_sizes);
                // timer.end();
                // timer.print("Write");                
            }

            write_metadata();
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
                            + get_size(level_elements)
                            + get_size(stopping_indices) + get_size(level_num) + get_size(decomposed_buffer_dims) 
                            + get_size(coeff_interp_directions) + get_size(structures) + 1; // one byte for whether negabinary encoding is used 
            uint8_t * metadata = (uint8_t *) malloc(metadata_size);
            uint8_t * metadata_pos = metadata;
            *(metadata_pos ++) = (uint8_t) dimensions.size();
            serialize(dimensions, metadata_pos);
            *(metadata_pos ++) = (uint8_t) level_error_bounds.size();
            serialize(level_error_bounds, metadata_pos);
            // serialize(level_squared_errors, metadata_pos);
            serialize(level_sizes, metadata_pos);
            serialize(level_elements, metadata_pos);
            serialize(stopping_indices, metadata_pos);
            serialize(level_num, metadata_pos);
            serialize(decomposed_buffer_dims, metadata_pos);
            serialize(coeff_interp_directions, metadata_pos);
            serialize(structures, metadata_pos);
            *(metadata_pos ++) = (uint8_t) negabinary;
            writer.write_metadata(metadata, metadata_size);
            free(metadata);
        }

        ~CPRefactor(){}

        void print() const {
            std::cout << "Coefficient prediction refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        bool refactor(uint8_t target_level, uint8_t num_bitplanes){
            uint8_t max_level = log2(*min_element(dimensions.begin(), dimensions.end())) - 1;
            if(target_level > max_level){
                std::cerr << "Target level is higher than " << max_level << std::endl;
                return false;
            }
            // Timer timer;
            // decompose data hierarchically
            // timer.start();
            T value_range = compute_value_range(data);
            std::vector<std::vector<T>> decomposed_buffers;
            if(target_level == 1){
                decomposed_buffers = decomposer.decompose_interleave(data.data(), dimensions, target_level);               
                decomposed_buffer_dims = decomposer.get_level_buffer_dims();    
            }
            else if (target_level > 1){
                auto decomposed_buffers_level_1 = decomposer.decompose_interleave(data.data(), dimensions, 1);
                decomposed_buffer_dims = decomposer.get_level_buffer_dims();
                // std::cout << decomposed_buffer_dims[0][0] << " " << decomposed_buffer_dims[0][1] << " " << decomposed_buffer_dims[0][2] << std::endl;
                decomposed_buffers = decomposer.decompose_interleave_combine_levels(decomposed_buffers_level_1[0].data(), decomposed_buffer_dims[0], target_level - 1);
                // std::cout << "decomposed_buffers.size() = " << decomposed_buffers.size() << "\n";
                for(int i=1; i<decomposed_buffers_level_1.size(); i++){
                    decomposed_buffers.push_back(decomposed_buffers_level_1[i]);
                }
                // std::cout << "decomposed_buffers.size() = " << decomposed_buffers.size() << "\n";
                decomposed_buffers_level_1.clear();
            }
            // auto decomposed_buffers = decomposer.decompose_interleave(data.data(), dimensions, target_level);
            // decomposed_buffer_dims = decomposer.get_level_buffer_dims();
            structures.resize(decomposed_buffer_dims.size() - 1);
            for(int i=0; i<structures.size(); i++){
                // structures[i].push_back(0);
                for(int j=0; j<target_level; j++){
                    structures[i].push_back(j);
                }
            }
            auto coeff_decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
            uint8_t coeff_target_level = 2;
            uint32_t coeff_stride = 15;
            uint32_t coeff_block_size = 9;

            // encode level by level
            level_error_bounds.clear();
            level_squared_errors.clear();
            level_components.clear();
            level_sizes.clear();
            std::vector<std::vector<T>> level_buffers;
            // level_buffers.push_back(decomposed_buffers[0]);
            for(int i=0; i<target_level; i++){
                level_buffers.push_back(decomposed_buffers[i]);
            }
            // // Encoder level 0 first
            // {
            //     int i=0;
            //     level_elements.push_back(level_buffers[i].size());
            //     T level_max_error = compute_max_abs_value(level_buffers[i].data(), level_buffers[i].size());
            //     if(negabinary) level_error_bounds.push_back(level_max_error * 4);
            //     else level_error_bounds.push_back(level_max_error);
            //     int level_exp = 0;
            //     frexp(level_max_error, &level_exp);
            //     std::vector<uint32_t> stream_sizes;
            //     auto streams = encoder.encode(level_buffers[i].data(), level_buffers[i].size(), level_exp, num_bitplanes, stream_sizes);
            //     uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
            //     stopping_indices.push_back(stopping_index);
            //     level_components.push_back(streams);
            //     level_sizes.push_back(stream_sizes);
            // }

            // Tune
            // std::cout << "Tuning" << std::endl;
            uint8_t tmp_num_level = target_level;
            for(int i=target_level; i<decomposed_buffers.size(); i++){
                // auto tuner = MDR::CoeffProfilingSamplingTuner<T, MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>, 
                //                                  Encoder, Compressor, MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, 
                //                                  MDR::MaxErrorEstimatorHB<T>>(coeff_decomposer, encoder, compressor, interpreter);
                // // tuner.copy_in_level_0_info(value_range, level_error_bounds[0], level_sizes[0]);
                // tuner.tune(decomposed_buffers[i].data(), decomposed_buffer_dims[i-target_level+1], coeff_target_level, num_bitplanes, coeff_stride, coeff_block_size);
                // coeff_interp_directions.push_back(tuner.get_best_direction());
                coeff_interp_directions.push_back(2);
                if(coeff_interp_directions.back() == -1) {
                    structures[i-target_level].push_back(tmp_num_level);
                    tmp_num_level++;
                }
                else {
                    for(int j=0; j<coeff_target_level+1; j++){
                        structures[i-target_level].push_back(tmp_num_level);
                        tmp_num_level++;
                    }
                }
            }
            for(int i=0; i<coeff_interp_directions.size(); i++){
                std::cout << "coeff_interp_directions[" << i << "] = " << (int)coeff_interp_directions[i] << std::endl;
            }

            // Coefficient Decomposition
            // std::cout << "Coefficient Decomposition" << std::endl;
            for(int i=target_level; i<decomposed_buffers.size(); i++){
                if(coeff_interp_directions[i-target_level] == -1){
                    level_buffers.push_back(decomposed_buffers[i]);
                }
                else{
                    coeff_decomposer.direction = coeff_interp_directions[i-target_level];
                    std::cout << "direction = " << (int)coeff_interp_directions[i - target_level] << std::endl;
                    std::cout << decomposed_buffer_dims[i-target_level+1][0] << " " << decomposed_buffer_dims[i-target_level+1][1] << " " << decomposed_buffer_dims[i-target_level+1][2] << std::endl;
                    // std::cout << "decomposed_buffers[" << i << "].size() = " << decomposed_buffers[i].size() << std::endl;
                    auto decomposed_coeff_buffers = coeff_decomposer.decompose_interleave_combine_levels(decomposed_buffers[i].data(), decomposed_buffer_dims[i-target_level+1], coeff_target_level);
                    for(int j=0; j<coeff_target_level+1; j++){
                        level_buffers.push_back(decomposed_coeff_buffers[j]);
                        // std::cout << "level_buffers[" << level_buffers.size() - 1 << "] = decomposed_coeff_buffers[" << j << "]" << std::endl;
                    }
                }
            }
            
            decomposed_buffers.clear();
            
            // std::cout << "decomposed data:" << std::endl;
            // for(int i=0; i<data.size(); i++){
            //     std::cout << data[i] << " ";
            // }
            // std::cout << std::endl;
            // MGARD::writefile("decomposed_coeff.dat", data.data(), data.size());
            // timer.end();
            // timer.print("Decompose");

            // auto level_dims = compute_level_dims_new(dimensions, target_level);
            // auto level_elements = compute_level_elements(level_dims, target_level);
            // std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
            // SquaredErrorCollector<T> s_collector = SquaredErrorCollector<T>();
            // std::cout << "Encoding" << std::endl;
            for(int i=0; i<level_buffers.size(); i++){
                // timer.start();
                // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                // T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
                // extract level i component
                // interleaver.interleave(data.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer), i, target_level);

                // std::cout << "Level " << i << " coefficients:" << std::endl;
                // for(int j=0; j < level_elements[i]; j++){
                //     std::cout << buffer[j] << " ";
                // }
                // std::cout << std::endl;
                // compute max coefficient as level error bound
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
                level_elements.push_back(level_buffers[i].size());
                T level_max_error = compute_max_abs_value(level_buffers[i].data(), level_buffers[i].size());
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
                std::vector<uint32_t> stream_sizes;
                // std::vector<double> level_sq_err;
                auto streams = encoder.encode(level_buffers[i].data(), level_buffers[i].size(), level_exp, num_bitplanes, stream_sizes);
                // free(buffer);
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
            std::cout << "level_error_bounds: " << std::endl;
            for(int i=0; i<level_error_bounds.size(); i++){
                std::cout << level_error_bounds[i] << " ";
            }
            std::cout << std::endl;
            return true;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        Compressor compressor;
        ErrorCollector collector;
        Writer writer;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<int> coeff_interp_directions;
        std::vector<std::vector<uint32_t>> decomposed_buffer_dims;
        std::vector<std::vector<uint8_t>> structures;
        std::vector<uint32_t> level_elements;
    public:
        bool negabinary = false;
    };
}
#endif

