#ifndef _MDR_PROFILING_SAMPLING_TUNER_HPP
#define _MDR_PROFILING_SAMPLING_TUNER_HPP

#include "TunerInterface.hpp"
#include "MDR/Decomposer/Decomposer.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/ErrorEstimator/ErrorEstimator.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/SizeInterpreter/SizeInterpreter.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/RefactorUtils.hpp"
#include "sample.hpp"
#include "qoi_utils.hpp"

namespace MDR{
    template<class T, class Level_Decomposer, class Layer_Decomposer, class Encoder, class Compressor, class Level_SizeInterpreter, class Layer_SizeInterpreter, class Level_ErrorEstimator, class Layer_ErrorEstimator>
    class ProfilingSamplingTuner : public concepts::TunerInterface<T> {
    public:
        ProfilingSamplingTuner(Level_Decomposer level_decomposer, Layer_Decomposer layer_decomposer, Encoder encoder, Compressor compressor, Level_SizeInterpreter level_interpreter, Layer_SizeInterpreter layer_interpreter)
            : level_decomposer(level_decomposer), layer_decomposer(layer_decomposer), encoder(encoder), compressor(compressor), level_interpreter(level_interpreter), layer_interpreter(layer_interpreter) {}
        void tune(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes, uint32_t stride=15, uint32_t block_size=5) {
            // std::cout << "Input: stride = " << (int) stride << ", block_size = " << (int) block_size << std::endl;
            if(block_size % 2 == 0){
                std::cout << "Even block_size = " << block_size << " is buggy due to the interpolation strategy, update block_size to " << 5 << std::endl;
                block_size = 5;
            }
            uint32_t num_elements = 1;
            for(int i=0; i<dims.size(); i++){
                num_elements *= dims[i];
                dimensions.push_back((size_t)dims[i]);
                // if(dims[i] / stride < 2){
                //     std::cout << "Dimension " << i << " = " << dims[i] << " is too small for sampling with stride = " << stride << ", update stride to " << dims[i] / 2 << std::endl; 
                //     stride = dims[i] / 2;
                // }
            }
            if(stride < block_size) {
                std::cout << "Sampling stride has to be greater than block_size" << std::endl;
                exit(-1);
            }
            // Timer timer;
            // timer.start();
            block_size = 4;
            for(int i=0; i<target_level; i++){
                block_size *= 2;
            }
            block_size --;
            uint32_t min_dim = *std::min_element(dims.begin(), dims.end());
            if(block_size > min_dim) block_size = min_dim - 1;
            // std::cout << "block_size = " << block_size << std::endl;
            std::vector<std::vector<size_t>> starts;
            MGARD::profiling_blocks<T>(data_, dimensions, starts, block_size, 1e-5, 1);
            MGARD::sample_blocks_after_profiling<T>(data_, dimensions, sampled_blocks, starts, block_size, 0.01);
            // std::cout << "one sampled_block size = " << sampled_blocks[0].size() << std::endl;
            // MGARD::sample_blocks<T>(data_, dimensions, sampled_blocks, (size_t)stride, (size_t)block_size);
            // std::cout << "sampled_blocks.size() = " << sampled_blocks.size() << std::endl;
            std::vector<double> ebs = {1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5};
            // std::vector<double> ebs = {1e-3};
            std::vector<std::vector<uint32_t>> retrieved_sizes(7);
            T value_range = compute_value_range(data_, num_elements);
            for(int i=0; i<ebs.size(); i++){
                ebs[i] *= value_range;
            }
            if(mode == 0){
                test_per_level_interp(target_level, num_bitplanes, block_size);
                // uint32_t total_retrieved_size = 0;
                for(auto tolerance : ebs){
                    uint32_t retrieved_size = test_per_level_reconstruct(tolerance, target_level + 1);
                    // std::cout << "Tolerance: " << tolerance << ", retrieved_size = " << retrieved_size << std::endl;
                    // total_retrieved_size += retrieved_size;
                    retrieved_sizes[0].push_back(retrieved_size);
                }
            }
            // std::cout << "Direction -1, total retrieved size = " << total_retrieved_size << std::endl;
            // uint32_t min_total_retrieved_size = total_retrieved_size;
            std::vector<std::vector<uint32_t>> interp_orders = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            int num_layers = target_level * dims.size() + 1;
            for(int i=0; i<interp_orders.size(); i++){
                layer_interpreter.reset();
                test_per_layer_interp(interp_orders[i], target_level, num_bitplanes, block_size);
                // uint32_t total_retrieved_size = 0;
                for(auto tolerance : ebs){
                    uint32_t retrieved_size = test_per_layer_reconstruct(tolerance, num_layers);
                    // std::cout << "Tolerance: " << tolerance << ", retrieved_size = " << retrieved_size << std::endl;
                    // total_retrieved_size += retrieved_size;
                    retrieved_sizes[i + 1].push_back(retrieved_size);
                }
                // std::cout << "Direction " << direction << ", total retrieved size = " << total_retrieved_size << std::endl;
                // if(total_retrieved_size < min_total_retrieved_size){
                //     min_total_retrieved_size = total_retrieved_size;
                //     best_direction = direction;
                // }
            }
            if(mode == 0){
                std::vector<int> votes(interp_orders.size() + 1, 0);
                for(int eb = 0; eb < ebs.size(); eb++){
                    uint32_t best_size = retrieved_sizes[0][eb];
                    int best_direct = 0;
                    for(int d = 1; d < interp_orders.size() + 1; d++){
                        if(retrieved_sizes[d][eb] < best_size){
                            best_size = retrieved_sizes[d][eb];
                            best_direct = d;
                        }
                    }
                    // std::cout << "best_direct = " << best_direct - 1 << std::endl;
                    votes[best_direct]++;
                }
                int best_interp = (std::max_element(votes.begin(), votes.end()) - votes.begin()) - 1;
                if(best_interp != -1){
                    best_interp_order = interp_orders[best_interp];
                }
            } else {
                std::vector<int> votes(interp_orders.size() + 1, 0);
                for(int eb = 0; eb < ebs.size(); eb++){
                    uint32_t best_size = std::numeric_limits<uint32_t>::max();
                    int best_direct = 0;
                    for(int d = 1; d < interp_orders.size() + 1; d++){
                        if(retrieved_sizes[d][eb] < best_size){
                            best_size = retrieved_sizes[d][eb];
                            best_direct = d;
                        }
                    }
                    // std::cout << "best_direct = " << best_direct - 1 << std::endl;
                    votes[best_direct]++;
                }
                int best_interp = (std::max_element(votes.begin(), votes.end()) - votes.begin()) - 1;
                // std::cout << "best_interp = " << best_interp << std::endl;
                best_interp_order = interp_orders[best_interp];
                // std::cout << "best_interp_order :" << std::endl;
                // for(int i=0; i<3; i++){
                //     std::cout << best_interp_order[i] << " ";
                // }
                // std::cout << std::endl;
            }
        }

        std::vector<uint32_t> get_best_interp_order(){
            return best_interp_order;
        }

        ~ProfilingSamplingTuner(){}

        void print() const {
            std::cout << "Profiling Sampling Tuner with the following components." << std::endl;
            std::cout << "Decomposer: "; level_decomposer.print(); layer_decomposer.print();
            std::cout << "Encoder:"; encoder.print();
            std::cout << "Interperter:"; level_interpreter.print(); layer_interpreter.print();
        }
    private:
        void test_per_level_interp(uint8_t target_level, uint8_t num_bitplanes, uint32_t block_size=5){
            level_error_bounds.clear();
            // level_error_bounds.push_back(level_0_error_bound);
            level_sizes.clear();
            // level_sizes.push_back(level_0_sizes);
            size_t num_dims = dimensions.size();
            std::vector<uint32_t> block_dims;
            std::vector<uint32_t> block_dims_uint32;
            for(int i=0; i<num_dims; i++){
                block_dims.push_back(block_size);
                block_dims_uint32.push_back(block_size);
            }
            auto level_dims = compute_level_dims_new(block_dims_uint32, target_level);
            auto level_elements = compute_level_elements(level_dims, target_level);
            // for(int i=0; i<level_dims.size(); i++){
            //     std::cout << "level " << i << ", dimensions = ";
            //     for(int j=0; j<level_dims[i].size(); j++){
            //         std::cout << level_dims[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // // for(int i=0; i<level_elements.size(); i++){
            //     std::cout << "level " << i << ", num of elements = " << level_elements[i] << std::endl;
            // }
            std::vector<std::vector<T>> level_buffers;
            std::vector<T*> level_buffer_ptr(target_level + 1);
            for(int i=0; i<=target_level; i++){
                level_buffers.push_back(std::vector<T>(level_elements[i] * sampled_blocks.size()));
                level_buffer_ptr[i] = level_buffers[i].data();
            }
            std::vector<std::vector<T>> sampled_blocks_copied = sampled_blocks;
            // decompose each block independently
            for(int b=0; b<sampled_blocks.size(); b++){
                // std::cout << compute_max_abs_value(sampled_blocks[b].data(), sampled_blocks[b].size()) << std::endl;
                auto block_level_buffers = level_decomposer.decompose_interleave_combine_levels(sampled_blocks_copied[b].data(), block_dims, target_level);
                // if(b == 0) std::cout << "block_level_buffers.size() = " << block_level_buffers.size() << std::endl;
                for(int i=0; i<=target_level; i++){
                    // if(b == 0) std::cout << block_level_buffers[i].size() << " " << level_elements[i] << std::endl;
                    memcpy(level_buffer_ptr[i], block_level_buffers[i].data(), level_elements[i] * sizeof(T));
                    level_buffer_ptr[i] += level_elements[i];
                }
            }
            // encoding together
            for(int i=0; i<= target_level; i++){
                T level_max_error = compute_max_abs_value(level_buffers[i].data(), level_buffers[i].size());
                if(negabinary) level_error_bounds.push_back(level_max_error * 4);
                else level_error_bounds.push_back(level_max_error);
                // std::cout << "level_max_error = " << level_max_error << std::endl;
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(level_buffers[i].data(), level_buffers[i].size(), level_exp, num_bitplanes, stream_sizes);
                uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                level_sizes.push_back(stream_sizes);
                // std::cout << "stream_sizes:" << std::endl;
                // for(int j=0; j<stream_sizes.size(); j++){
                //     std::cout << (int) stream_sizes[i] << " ";
                // }
                // std::cout << std::endl;
            }
        }
        void test_per_layer_interp(std::vector<uint32_t> interp_order, uint8_t target_level, uint8_t num_bitplanes, uint32_t block_size=5){
            level_error_bounds.clear();
            // level_error_bounds.push_back(level_0_error_bound);
            level_sizes.clear();
            // level_sizes.push_back(level_0_sizes);
            layer_decomposer.interp_order = interp_order;
            size_t num_dims = dimensions.size();
            std::vector<uint32_t> block_dims;
            std::vector<uint32_t> block_dims_uint32;
            for(int i=0; i<num_dims; i++){
                block_dims.push_back(block_size);
                block_dims_uint32.push_back(block_size);
            }
            auto level_dims = MGARD::compute_level_dims_new(block_dims_uint32, target_level);
            std::vector<std::vector<uint32_t>> level_buffer_dims;
            auto level_elements = MGARD::compute_level_buffers_size_generic(level_dims, target_level, interp_order, level_buffer_dims);
            // for(int i=0; i<level_dims.size(); i++){
            //     std::cout << "level " << i << ", dimensions = ";
            //     for(int j=0; j<level_dims[i].size(); j++){
            //         std::cout << level_dims[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // for(int i=0; i<level_elements.size(); i++){
            //     std::cout << "level " << i << ", num of elements = " << level_elements[i] << std::endl;
            // }
            std::vector<std::vector<T>> level_buffers;
            int num_layers = target_level * num_dims + 1;
            std::vector<T*> level_buffer_ptr(num_layers);
            for(int i=0; i<num_layers; i++){
                level_buffers.push_back(std::vector<T>(level_elements[i] * sampled_blocks.size()));
                level_buffer_ptr[i] = level_buffers[i].data();
                // std::cout << level_buffers[i].size() << std::endl;
            }
            std::vector<std::vector<T>> sampled_blocks_copied = sampled_blocks;
            // decompose each block independently
            for(int b=0; b<sampled_blocks.size(); b++){
                // std::cout << compute_max_abs_value(sampled_blocks[b].data(), sampled_blocks[b].size()) << std::endl;
                auto block_level_buffers = layer_decomposer.decompose_interleave(sampled_blocks_copied[b].data(), block_dims, target_level);
                // if(b == 0) std::cout << "block_level_buffers.size() = " << block_level_buffers.size() << std::endl;
                for(int i=0; i<num_layers; i++){
                    // if(b == 0) std::cout << block_level_buffers[i].size() << " " << level_elements[i] << std::endl;
                    memcpy(level_buffer_ptr[i], block_level_buffers[i].data(), level_elements[i] * sizeof(T));
                    level_buffer_ptr[i] += level_elements[i];
                }
            }
            // encoding together
            for(int i=0; i<num_layers; i++){
                T level_max_error = compute_max_abs_value(level_buffers[i].data(), level_buffers[i].size());
                if(negabinary) level_error_bounds.push_back(level_max_error * 4);
                else level_error_bounds.push_back(level_max_error);
                // std::cout << "level_max_error = " << level_max_error << std::endl;
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(level_buffers[i].data(), level_buffers[i].size(), level_exp, num_bitplanes, stream_sizes);
                uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                level_sizes.push_back(stream_sizes);
                // std::cout << "stream_sizes:" << std::endl;
                // for(int j=0; j<stream_sizes.size(); j++){
                //     std::cout << (int) stream_sizes[i] << " ";
                // }
                // std::cout << std::endl;
            }
        }

        uint32_t test_per_level_reconstruct(double tolerance, uint8_t num_levels){
            std::vector<std::vector<double>> level_errors;
            MaxErrorCollector<T> collector = MaxErrorCollector<T>();
            for(int i=0; i<num_levels; i++){
                auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                level_errors.push_back(collected_error);
            }
            std::vector<uint8_t> level_num_bitplanes(num_levels, 0);
            // std::cout << "level_error_bounds: " << std::endl;
            // for(int i=0; i<= target_level; i++){
            //     std::cout << level_error_bounds[i] << " ";
            // }
            // std::cout << std::endl;
            auto retrieve_sizes = level_interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
            // std::cout << "level_num_bitplanes: " << std::endl;
            // for(int i=0; i<= target_level; i++){
            //     std::cout << (int) level_num_bitplanes[i] << " ";
            // }
            // std::cout << std::endl;
            uint32_t retrieved_size = 0;
            for(int i=0; i<retrieve_sizes.size(); i++){
                retrieved_size += retrieve_sizes[i];
            }
            return retrieved_size;
        }

        uint32_t test_per_layer_reconstruct(double tolerance, uint8_t num_levels){
            std::vector<std::vector<double>> level_errors;
            MaxErrorCollector<T> collector = MaxErrorCollector<T>();
            for(int i=0; i<num_levels; i++){
                auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                level_errors.push_back(collected_error);
            }
            std::vector<uint8_t> level_num_bitplanes(num_levels, 0);
            // std::cout << "level_error_bounds: " << std::endl;
            // for(int i=0; i<= target_level; i++){
            //     std::cout << level_error_bounds[i] << " ";
            // }
            // std::cout << std::endl;
            auto retrieve_sizes = layer_interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
            // std::cout << "level_num_bitplanes: " << std::endl;
            // for(int i=0; i<= target_level; i++){
            //     std::cout << (int) level_num_bitplanes[i] << " ";
            // }
            // std::cout << std::endl;
            uint32_t retrieved_size = 0;
            for(int i=0; i<retrieve_sizes.size(); i++){
                retrieved_size += retrieve_sizes[i];
            }
            return retrieved_size;
        }

        Level_Decomposer level_decomposer;
        Layer_Decomposer layer_decomposer;
        Encoder encoder;
        Compressor compressor;
        Level_SizeInterpreter level_interpreter;
        Layer_SizeInterpreter layer_interpreter;
        std::vector<std::vector<T>> sampled_blocks;
        std::vector<size_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> best_interp_order = {0, 0, 0};
        T value_range;
        T level_0_error_bound;
        std::vector<uint32_t> level_0_sizes;
    public:
        uint8_t mode = 0; // bitrate prior = 0, psnr prior = 1
        bool negabinary = false;
    };
}
#endif