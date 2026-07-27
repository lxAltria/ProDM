#ifndef _MDR_ORDERED_COEFFICIENT_PREDICTION_REFACTOR_2_HPP
#define _MDR_ORDERED_COEFFICIENT_PREDICTION_REFACTOR_2_HPP

#include "RefactorInterface.hpp"
#include "MDR/Decomposer/Decomposer.hpp"
#include "MDR/Interleaver/Interleaver.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/ErrorEstimator/ErrorEstimator.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/Writer/Writer.hpp"
#include "MDR/RefactorUtils.hpp"
#include "MDR/Tuner/Tuner.hpp"
#include "MDR/SizeInterpreter/SizeInterpreter.hpp"
#include <queue>

namespace MDR {
    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, and error collector
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
    class OrderedCPRefactor_2 : public concepts::RefactorInterface<T> {
    public:
        OrderedCPRefactor_2(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), collector(collector), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){
            Timer timer;
            timer.start();
            dimensions = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            double value_range = compute_value_range(data);
            // if refactor successfully
            if(refactor(target_level, num_bitplanes)){
                timer.end();
                timer.print("Refactor");
                // timer.start();
                // level_num = writer.write_level_components(level_components, level_sizes);
                // timer.end();
                // timer.print("Write");                
            }

            // Getting error table
            std::vector<std::vector<double>> level_abs_errors;
            std::vector<std::vector<double>>& level_errors = level_squared_errors;
            {
                std::cout << "Computing absolute error" << std::endl;
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for(int i=0; i<level_sizes.size(); i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }

            std::vector<std::vector<std::vector<uint8_t>>> layer_order(coeff_interp_levels, std::vector<std::vector<uint8_t>>());
            std::vector<std::vector<std::vector<double>>> layer_error_perstep(coeff_interp_levels, std::vector<std::vector<double>>());
            std::vector<std::vector<double>> combined_coeff_error_perstep(coeff_interp_levels, std::vector<double>());
            coeff_components.resize(coeff_interp_levels);
            coeff_sizes.resize(coeff_interp_levels);
            coeff_combined_bps.resize(coeff_interp_levels);
            uint8_t start_level = target_level - coeff_interp_levels + 1;
            double max_region_error = 0;
            auto compute_region_max_error = [&](const std::vector<std::vector<double>>& level_errors_, uint8_t start_level, uint8_t num_levels){
                double sum = 0;
                for(uint8_t i = start_level; i < start_level + num_levels; i++){
                    sum += level_errors_[i][0];
                }
                return sum;
            };

            if(greedy_bfs){
                std::cout << "Greedy + max + BFS" << std::endl;
                for(int i=0; i<coeff_interp_levels; i++){
                    for(int j=0; j<dims.size(); j++){
                        // std::cout << "start_level = " << (int) start_level << std::endl;
                        int region_idx = i * dims.size() + j;
                        std::vector<double> tmp_error_perstep;
                        std::vector<uint8_t> tmp_order;
                        if(coeff_interp_directions[region_idx] == -1){
                            if(max_region_error < level_errors[start_level][0]) max_region_error = level_errors[start_level][0];
                            tmp_order = get_chunks_order_greedy_coefficient(level_errors, start_level, 1, tmp_error_perstep);
                            tmp_error_perstep.insert(tmp_error_perstep.begin(), level_errors[start_level][0]);
                            start_level++;
                        } 
                        else {
                            double tmp_max_region_error = compute_region_max_error(level_errors, start_level, coeff_target_level + 1);
                            if(max_region_error < tmp_max_region_error) max_region_error = tmp_max_region_error;
                            tmp_order = get_chunks_order_greedy_coefficient(level_errors, start_level, coeff_target_level + 1, tmp_error_perstep);
                            tmp_error_perstep.insert(tmp_error_perstep.begin(), tmp_max_region_error);
                            start_level += (coeff_target_level + 1);
                        }
                        layer_order[i].push_back(tmp_order);
                        layer_error_perstep[i].push_back(tmp_error_perstep);
                    }
                    std::vector<double> coefficient_error_perstep;
                    std::vector<uint8_t> coefficient_order = get_chunks_order_max_coeffcient(layer_error_perstep[i], layer_order[i], coefficient_error_perstep);
                    combined_coeff_error_perstep[i] = combine_coeff_bitplanes_after_max_order(coefficient_error_perstep, coefficient_order, i);
                    double error_propagate = 1;
                    for(int j=0; j<coeff_interp_levels-i-1; j++){
                        error_propagate = error_propagate * 1.25 * 1.25 * 1.25;
                    }
                    for(int j=0; j<combined_coeff_error_perstep[i].size(); j++){
                        combined_coeff_error_perstep[i][j] *= error_propagate;
                    }
                    // std::cout << "combined_coeff_error_perstep[" << i << "].size() = " << combined_coeff_error_perstep[i].size() << std::endl;
                }
                // combined_coeff_error_perstep.insert(combined_coeff_error_perstep.begin(), max_region_error);
                auto overall_estimator = MDR::MaxErrorEstimatorHBCubic<T>(dims.size());
                auto overall_interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(overall_estimator);
                std::vector<double> tmp_tolerances = {5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10, 5e-11, 1e-11, 0};
                // std::vector<double> tmp_tolerances = {5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 0};
                std::vector<uint8_t> tmp_index(target_level - coeff_interp_levels + 1, 0);
                std::vector<uint16_t> coeff_index(coeff_interp_levels);
                for(int i=0; i<tmp_tolerances.size(); i++){
                    tmp_tolerances[i] *= value_range;
                    auto output = overall_interpreter.interpret_overall_retrieve_size_2(level_sizes, level_errors, coeff_sizes, combined_coeff_error_perstep, target_level - coeff_interp_levels + 1, tmp_tolerances[i], tmp_index, coeff_index);
                }
                level_chunk_order = overall_interpreter.get_path();
                level_error_perstep = overall_interpreter.get_error_perstep();
            }

            write_metadata();

            // Calculate total size
            start_level = target_level - coeff_interp_levels + 1;
            std::cout << std::endl;
            std::vector<int> consumed(target_level + 1, 0);
            uint64_t total_size_64 = 0;
            std::vector<size_t> idx_coeff_combined_bps(coeff_interp_levels, 0);
            for (uint8_t lev : level_chunk_order) {
                int j = consumed[lev];
                if(lev < start_level) total_size_64 += level_sizes[lev][j];
                else {
                    int coeff_l = lev - start_level;
                    total_size_64 += coeff_sizes[coeff_l][idx_coeff_combined_bps[coeff_l]];
                    idx_coeff_combined_bps[coeff_l]++;
                }
                consumed[lev] = j + 1;
            }

            // Attention: file size cannot be over 4GB
            uint32_t total_size = static_cast<uint32_t>(total_size_64);

            std::vector<uint8_t> packed(total_size);
            uint8_t* dst = packed.data();
            // copy every chunk
            std::fill(consumed.begin(), consumed.end(), 0);
            std::fill(idx_coeff_combined_bps.begin(), idx_coeff_combined_bps.end(), 0);
            for (uint8_t lev : level_chunk_order) {
                int j = consumed[lev];
                if(lev < start_level){
                    uint32_t sz = level_sizes[lev][j];
                    std::memcpy(dst, level_components[lev][j], sz);
                    dst += sz;
                } else {
                    int coeff_l = lev - start_level;
                    uint32_t sz = coeff_sizes[coeff_l][idx_coeff_combined_bps[coeff_l]];
                    std::memcpy(dst, coeff_components[coeff_l][idx_coeff_combined_bps[coeff_l]].data(), sz);
                    idx_coeff_combined_bps[coeff_l]++;
                    dst += sz;
                }
                consumed[lev] = j + 1;
            }

            writer.write_components(packed.data(), total_size);

            for(int i=0; i<level_components.size(); i++){
                for(int j=0; j<level_components[i].size(); j++){
                    free(level_components[i][j]);                    
                }
            }
        }

        uint8_t * get_metadata(uint32_t& metadata_size) const {
            metadata_size =
                sizeof(uint8_t)  + get_size(dimensions)
                + sizeof(uint8_t)  + get_size(level_error_bounds)  
                + sizeof(uint8_t)   // interpolation type: 0 - per level, 1 - per layer
                + get_size(level_sizes)                            
                + get_size(level_elements)
                + get_size(stopping_indices)
                + sizeof(uint8_t)    // negabinary
                + sizeof(uint16_t)   // chunk_num
                + get_size(level_chunk_order)                                     
                + get_size(level_error_perstep)
                + sizeof(uint8_t)// coeff_interp_levels
                + get_size(decomposed_buffer_dims)
                + get_size(coeff_interp_directions); 
            for(int i=0; i<coeff_interp_levels; i++){
                metadata_size += sizeof(uint16_t) + get_size(coeff_combined_bps[i]);
            }
                
            uint8_t* metadata = static_cast<uint8_t*>(malloc(metadata_size));
            uint8_t* p = metadata;

            *(p++) = (uint8_t) dimensions.size();
            serialize(dimensions, p);
            *(p++) = (uint8_t) level_error_bounds.size();
            serialize(level_error_bounds, p);
            *(p++) = static_cast<uint8_t>(0);
            serialize(level_sizes, p);
            serialize(level_elements, p);
            serialize(stopping_indices, p);
            *(p++) = static_cast<uint8_t>(negabinary);

            const uint16_t chunk_num = level_chunk_order.size();
            memcpy(p, &chunk_num, sizeof(uint16_t));
            p += sizeof(uint16_t);
            serialize(level_chunk_order, p);
            serialize(level_error_perstep, p);
            *(p++) = static_cast<uint8_t>(coeff_interp_levels);
            serialize(decomposed_buffer_dims, p);
            serialize(coeff_interp_directions, p);
            for(int i=0; i<coeff_interp_levels; i++){
                const uint16_t coeff_combined_bp_num = coeff_combined_bps[i].size();
                memcpy(p, &coeff_combined_bp_num, sizeof(uint16_t));
                // std::cout << "coeff_combined_bp_num = " << (int)coeff_combined_bp_num << std::endl;
                p += sizeof(uint16_t);
                serialize(coeff_combined_bps[i], p);
            }

            return metadata;
        }

        void write_metadata() const {
            uint32_t metadata_size;
            uint8_t* metadata = get_metadata(metadata_size);

            writer.write_metadata(metadata, metadata_size);
            free(metadata);
        }

        ~OrderedCPRefactor_2(){}

        void print() const {
            std::cout << "Ordered coefficient prediction refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        inline bool approx_equal(double a, double b) const {
            constexpr double rel = 1e-3;
            double scale = std::max(std::abs(a), std::abs(b));
            return std::abs(a - b) <= std::min(rel * scale, 1e-12);
        }

        std::vector<uint8_t> get_chunks_order_greedy_coefficient(const std::vector<std::vector<double>>& level_errors, const uint8_t start_level, const uint8_t num_levels, std::vector<double>& error_perstep) const {
            // for(int i=start_level; i<start_level + num_levels; i++){
            //     for(int j=0; j<level_errors[i].size(); j++){
            //         std::cout << level_errors[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            auto error_estimator = MDR::MaxErrorEstimatorHB<T>();
            std::vector<uint8_t> index(level_sizes.size(), 0);
            int end_level = start_level + num_levels;
            double accumulated_error = 0;
            for(int i=start_level; i<end_level; i++){
                accumulated_error += error_estimator.estimate_error(level_errors[i][index[i]], i);
            }
            std::priority_queue<UnitErrorGain, std::vector<UnitErrorGain>, CompareUnitErrorGain> heap;
            // identify minimal level
            double min_error = accumulated_error;

            std::vector<uint8_t> order;
            std::vector<std::vector<uint32_t>> level_sizes_used(level_sizes.begin() + start_level, level_sizes.begin() + end_level);
            order.reserve([&](){
                size_t s=0; for(const auto& v: level_sizes_used) s+=v.size(); return s;
            }()); // allocate enough space for it

            for(int i=start_level; i<end_level; i++){
                min_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
                min_error += error_estimator.estimate_error(level_errors[i].back(), i);
                // fetch the first component if index is 0
                if(index[i] == 0){
                    accumulated_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
                    accumulated_error += error_estimator.estimate_error(level_errors[i][index[i] + 1], i);
                    index[i] ++;
                    order.push_back(static_cast<uint8_t>(i));
                    error_perstep.push_back(accumulated_error);
                    // std::cout << i << " ";
                }
                // push the next one
                if(index[i] != level_sizes[i].size()){
                    double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i] + 1], i);
                    heap.push(UnitErrorGain(error_gain / level_sizes[i][index[i]], i));
                }
            }

            while(!heap.empty()){
                auto unit_error_gain = heap.top();
                heap.pop();
                int i = unit_error_gain.level;
                int j = index[i];
                accumulated_error -= error_estimator.estimate_error(level_errors[i][j], i);
                accumulated_error += error_estimator.estimate_error(level_errors[i][j + 1], i);
                index[i] ++;
                if(index[i] != level_sizes[i].size()){
                    double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i] + 1], i);
                    heap.push(UnitErrorGain(error_gain / level_sizes[i][index[i]], i));
                }
                order.push_back(static_cast<uint8_t>(i));
                error_perstep.push_back(accumulated_error);
            }
            // std::cout << "Region error_perstep " << std::endl;
            // for(int i=0; i<error_perstep.size(); i++){
            //     std::cout << "(" << (int)order[i] << ", " << error_perstep[i] << ") ";
            // }
            // std::cout << std::endl;
            return order;
        }

        std::vector<uint8_t> get_chunks_order_max_coeffcient(const std::vector<std::vector<double>>& layer_error_perstep, const std::vector<std::vector<uint8_t>>& layer_order, std::vector<double>& error_perstep) const {
            // for(int i=0; i<layer_error_perstep.size(); i++){
            //     for(int j=0; j<layer_error_perstep[i].size(); j++){
            //         std::cout << layer_error_perstep[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            std::vector<uint8_t> order;
            order.resize([&]{
                size_t s = 0; for(const auto & v: layer_order) s+=v.size(); return s;
            }());
            error_perstep.resize(order.size());
            size_t count = 0;
            std::vector<size_t> index(layer_order.size(), 0);
            while(count < order.size()){
                int max_layer = -1;
                double max_error = std::numeric_limits<double>::lowest();
                for(int i=0; i<layer_order.size(); i++){
                    if(index[i] >= layer_order[i].size()) continue;
                    if(max_error < layer_error_perstep[i][index[i]]){
                        max_layer = i;
                        max_error = layer_error_perstep[i][index[i]];
                    }
                }
                order[count] = layer_order[max_layer][index[max_layer]];
                error_perstep[count] = layer_error_perstep[max_layer][index[max_layer]];
                count++;
                index[max_layer]++;
            }
            // for(int i=0; i<count; i++){
            //     std::cout << "(" << (int)order[i] << ", " << error_perstep[i] << ") ";
            // }
            // std::cout << std::endl;
            return order;
        }

        std::vector<double> combine_coeff_bitplanes_after_max_order(const std::vector<double>& error_perstep, const std::vector<uint8_t>& order, int level_coeff){
            std::vector<double> coeff_error_perstep;
            int idx = 0;
            std::vector<size_t> index(level_sizes.size(), 0);
            // get coeff_combined_bps & coeff_sizes & coeff_error_perstep
            while(idx < error_perstep.size()){
                coeff_error_perstep.push_back(error_perstep[idx]);
                std::vector<uint8_t> tmp_bps;
                if((idx + 1 < error_perstep.size()) && (approx_equal(error_perstep[idx], error_perstep[idx + 1]))){
                    tmp_bps.push_back(order[idx]);
                    uint32_t tmp_size = level_sizes[order[idx]][index[order[idx]]++];
                    while((idx + 1 < error_perstep.size()) && (approx_equal(error_perstep[idx], error_perstep[idx + 1]))){
                        idx++;
                        tmp_bps.push_back(order[idx]);
                        tmp_size += level_sizes[order[idx]][index[order[idx]]++];
                    }
                    coeff_sizes[level_coeff].push_back(tmp_size);
                } else {
                    tmp_bps.push_back(order[idx]);
                    coeff_sizes[level_coeff].push_back(level_sizes[order[idx]][index[order[idx]]++]);
                }
                coeff_combined_bps[level_coeff].push_back(tmp_bps);
                idx++;
            }
            // /*
            // size_t bps_count = 0;
            // for(int i=0; i<coeff_combined_bps.size(); i++){
            //     bps_count += coeff_combined_bps[i].size();
            // }
            // std::cout << "bps_count = " << bps_count << std::endl;
            // uint32_t coeff_total_size = 0;
            // for(int i=0; i<coeff_sizes.size(); i++){
            //     coeff_total_size += coeff_sizes[i];
            // }
            // std::cout << "coeff_total_size = " << coeff_total_size << std::endl;
            // uint32_t coeff_levels_total_size = 0;
            // for(int i=4; i<13; i++){
            //     for(int j=0; j<level_sizes[i].size(); j++){
            //         coeff_levels_total_size += level_sizes[i][j];
            //     }
            // }
            // std::cout << "coeff_levels_total_size = " << coeff_levels_total_size << std::endl;
            // int step = 0;
            // for(int i=0; i<coeff_error_perstep.size(); i++){
            //     std::cout << "(";
            //     for(int j=0; j<coeff_combined_bps[i].size(); j++){
            //         std::cout << (int) coeff_combined_bps[i][j] << ", ";
            //     }
            //     std::cout << coeff_error_perstep[i] << ", " << coeff_sizes[i] << ") ";
            //     step++;
            //     if(step == 5){
            //         step = 0;
            //         std::cout << std::endl;
            //     }
            // }
            // std::cout << std::endl;
            // */

            // malloc space for memcpy from level_components to coeff_components
            for(int i=0; i<coeff_sizes[level_coeff].size(); i++){
                coeff_components[level_coeff].push_back(std::vector<uint8_t>(coeff_sizes[level_coeff][i], 0));
            }
            std::fill(index.begin(), index.end(), 0);
            for(int i=0; i<coeff_combined_bps[level_coeff].size(); i++){
                uint8_t * dst = coeff_components[level_coeff][i].data();
                for(int j=0; j<coeff_combined_bps[level_coeff][i].size(); j++){
                    memcpy(dst, level_components[coeff_combined_bps[level_coeff][i][j]][index[coeff_combined_bps[level_coeff][i][j]]], level_sizes[coeff_combined_bps[level_coeff][i][j]][index[coeff_combined_bps[level_coeff][i][j]]]);
                    dst += level_sizes[coeff_combined_bps[level_coeff][i][j]][index[coeff_combined_bps[level_coeff][i][j]]++];
                }
            }
            // std::cout << "index" << std::endl;
            // for(int i=0; i<index.size(); i++){
            //     std::cout << (int)index[i] << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "coeff_combined_bps.size() = " << coeff_combined_bps.size() << ", coeff_sizes.size() = " << coeff_sizes.size() << ", coeff_error_perstep.size() = " << coeff_error_perstep.size() << ", coeff_components.size() = " << coeff_components.size() << std::endl;
            // std::cout << "From combine_coeff_bitplanes_after_max_order: coeff_error_perstep[" << 63 << "] = " << coeff_error_perstep[63] << std::endl;
            return coeff_error_perstep;
        }

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
            if(target_level == coeff_interp_levels){
                decomposed_buffers = decomposer.decompose_interleave(data.data(), dimensions, target_level);               
                decomposed_buffer_dims = decomposer.get_level_buffer_dims();    
            }
            else if (target_level > coeff_interp_levels){
                auto decomposed_buffers_level_1 = decomposer.decompose_interleave(data.data(), dimensions, coeff_interp_levels);
                decomposed_buffer_dims = decomposer.get_level_buffer_dims();
                std::cout << "decomposed_buffer_dims.size() = " << decomposed_buffer_dims.size() << std::endl;
                // std::cout << decomposed_buffer_dims[0][0] << " " << decomposed_buffer_dims[0][1] << " " << decomposed_buffer_dims[0][2] << std::endl;
                decomposed_buffers = decomposer.decompose_interleave_combine_levels(decomposed_buffers_level_1[0].data(), decomposed_buffer_dims[0], target_level - coeff_interp_levels);
                // std::cout << "decomposed_buffers.size() = " << decomposed_buffers.size() << "\n";
                for(int i=1; i<decomposed_buffers_level_1.size(); i++){
                    decomposed_buffers.push_back(decomposed_buffers_level_1[i]);
                }
                // std::cout << "decomposed_buffers.size() = " << decomposed_buffers.size() << "\n";
            }
            // auto decomposed_buffers = decomposer.decompose_interleave(data.data(), dimensions, target_level);
            // decomposed_buffer_dims = decomposer.get_level_buffer_dims();
            
            auto coeff_decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            auto coeff_interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
            uint32_t coeff_stride = 15;
            uint32_t coeff_block_size = 9;
            int no_cp_levels = target_level + 1 - coeff_interp_levels; // 4 + 1 - 2 = 3

            // encode level by level
            level_error_bounds.clear();
            level_squared_errors.clear();
            level_components.clear();
            level_sizes.clear();
            std::vector<std::vector<T>> level_buffers;
            // level_buffers.push_back(decomposed_buffers[0]);
            for(int i=0; i<no_cp_levels; i++){
                level_buffers.push_back(decomposed_buffers[i]);
            }

            // Tune
            if(!coeff_interp_directions.size()){
                std::cout << "Tuning" << std::endl;
                for(int i=no_cp_levels; i<decomposed_buffers.size(); i++){
                    auto tuner = MDR::CoeffProfilingSamplingTuner<T, MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>, 
                                                     Encoder, Compressor, MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, 
                                                     MDR::MaxErrorEstimatorHB<T>>(coeff_decomposer, encoder, compressor, coeff_interpreter);
                    tuner.tune(decomposed_buffers[i].data(), decomposed_buffer_dims[i-no_cp_levels+1], coeff_target_level, num_bitplanes, coeff_stride, coeff_block_size);
                    coeff_interp_directions.push_back(tuner.get_best_direction());
                    // coeff_interp_directions.push_back(2);
                }
            }
            std::cout << "coeff_interp_directions = ";
            for(int i=0; i<coeff_interp_directions.size(); i++){
                std::cout << (int)coeff_interp_directions[i] << " ";
            }
            std::cout << std::endl;

            // Coefficient Decomposition
            // std::cout << "Coefficient Decomposition" << std::endl;
            for(int i=no_cp_levels; i<decomposed_buffers.size(); i++){
                coeff_value_ranges.push_back(compute_value_range(decomposed_buffers[i]));
                if(coeff_interp_directions[i-no_cp_levels] == -1){
                    level_buffers.push_back(decomposed_buffers[i]);
                }
                else{
                    coeff_decomposer.direction = coeff_interp_directions[i-no_cp_levels];
                    // std::cout << "direction = " << (int)coeff_interp_directions[i - target_level] << std::endl;
                    std::cout << decomposed_buffer_dims[i-no_cp_levels+1][0] << " " << decomposed_buffer_dims[i-no_cp_levels+1][1] << " " << decomposed_buffer_dims[i-no_cp_levels+1][2] << std::endl;
                    // std::cout << "decomposed_buffers[" << i << "].size() = " << decomposed_buffers[i].size() << std::endl;
                    auto decomposed_coeff_buffers = coeff_decomposer.decompose_interleave_combine_levels(decomposed_buffers[i].data(), decomposed_buffer_dims[i-no_cp_levels+1], coeff_target_level);
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
            // std::cout << "level_sizes: " << std::endl;
            // for(int i=0; i<level_sizes.size(); i++){
            //     std::cout << "level " << i << std::endl;
            //     uint32_t level_x_total_size = 0;
            //     for(int j=0; j<level_sizes[i].size(); j++){
            //         level_x_total_size += level_sizes[i][j];
            //         std::cout << level_sizes[i][j] << " ";
            //     }
            //     std::cout << "\nTotal size = " << level_x_total_size << std::endl;
            // }
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
        std::vector<std::vector<uint32_t>> decomposed_buffer_dims;
        std::vector<uint32_t> level_elements;
        std::vector<uint8_t> level_chunk_order;
        std::vector<double> level_error_perstep;
        std::vector<std::vector<std::vector<uint8_t>>> coeff_components;
        std::vector<std::vector<uint32_t>> coeff_sizes;
        std::vector<std::vector<std::vector<uint8_t>>> coeff_combined_bps;
        uint8_t coeff_target_level = 2;
        std::vector<double> coeff_value_ranges;
    public:
        bool negabinary = false;
        bool greedy = false;
        bool bfs = false;
        bool greedy_bfs = true;
        std::vector<int> coeff_interp_directions;
        int coeff_interp_levels = 1;
    };
}
#endif

