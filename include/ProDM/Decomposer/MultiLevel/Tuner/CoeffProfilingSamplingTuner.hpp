#ifndef _MDR_COEFFICIENT_PROFILING_SAMPLING_TUNER_HPP
#define _MDR_COEFFICIENT_PROFILING_SAMPLING_TUNER_HPP

#include "TunerInterface.hpp"
#include "ProDM/Decomposer/MultiLevel/Decomposer.hpp"
#include "ProDM/Encoder/BitplaneEncoder.hpp"
#include "ProDM/ErrorControl/Estimator/ErrorEstimator.hpp"
#include "ProDM/ErrorControl/Collector/ErrorCollector.hpp"
#include "ProDM/ErrorControl/SizeInterpreter/SizeInterpreter.hpp"
#include "ProDM/Compressor/LevelCompressor.hpp"
#include "ProDM/Utils/RefactorUtils.hpp"
#include "ProDM/Decomposer/MultiLevel/MGARDx/sample.hpp"
#include "ProDM/Utils/QoIUtils.hpp"

namespace MDR{
    template<class T, class Decomposer, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator>
    class CoeffProfilingSamplingTuner : public concepts::TunerInterface<T> {
    public:
        CoeffProfilingSamplingTuner(Decomposer decomposer, Encoder encoder, Compressor compressor, SizeInterpreter interpreter)
            : decomposer(decomposer), encoder(encoder), compressor(compressor), interpreter(interpreter) {}
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
            std::vector<std::vector<size_t>> starts;
            MGARD::profiling_blocks<T>(data_, dimensions, starts, block_size, 1e-5, 1);
            MGARD::sample_blocks_after_profiling<T>(data_, dimensions, sampled_blocks, starts, block_size, 0.01);
            // MGARD::sample_blocks<T>(data_, dimensions, sampled_blocks, (size_t)stride, (size_t)block_size);
            // std::cout << "sampled_blocks.size() = " << sampled_blocks.size() << std::endl;
            // std::vector<double> ebs = {1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5};
            std::vector<double> ebs = {1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7};
            // std::vector<double> ebs = {1e-3};
            std::vector<std::vector<uint32_t>> retrieved_sizes(4);
            T value_range = compute_value_range(data_, num_elements);
            test_direct_BP_encoding(num_bitplanes, block_size);
            for(int i=0; i<ebs.size(); i++){
                ebs[i] *= value_range;
            }
            // uint32_t total_retrieved_size = 0;
            for(auto tolerance : ebs){
                uint32_t retrieved_size = test_reconstruct(tolerance, 0);
                // std::cout << "Tolerance: " << tolerance << ", retrieved_size = " << retrieved_size << std::endl;
                // total_retrieved_size += retrieved_size;
                retrieved_sizes[0].push_back(retrieved_size);
            }
            // std::cout << "Direction -1, total retrieved size = " << total_retrieved_size << std::endl;
            // uint32_t min_total_retrieved_size = total_retrieved_size;
            best_direction = -1;
            for(int direction = 0; direction < 3; direction++){
                test_refactor(direction, target_level, num_bitplanes, block_size);
                // uint32_t total_retrieved_size = 0;
                for(auto tolerance : ebs){
                    uint32_t retrieved_size = test_reconstruct(tolerance, target_level);
                    // std::cout << "Tolerance: " << tolerance << ", retrieved_size = " << retrieved_size << std::endl;
                    // total_retrieved_size += retrieved_size;
                    retrieved_sizes[direction + 1].push_back(retrieved_size);
                }
                // std::cout << "Direction " << direction << ", total retrieved size = " << total_retrieved_size << std::endl;
                // if(total_retrieved_size < min_total_retrieved_size){
                //     min_total_retrieved_size = total_retrieved_size;
                //     best_direction = direction;
                // }
            }
            std::vector<int> votes(4, 0);
            for(int eb = 0; eb < ebs.size(); eb++){
                uint32_t best_size = retrieved_sizes[0][eb];
                int best_direct = 0;
                for(int d = 1; d < 4; d++){
                    if(retrieved_sizes[d][eb] < best_size){
                        best_size = retrieved_sizes[d][eb];
                        best_direct = d;
                    }
                }
                // std::cout << "best_direct = " << best_direct - 1 << std::endl;
                votes[best_direct]++;
            }
            // if (votes[1] + votes[2] + votes[3] > 0){
            best_direction = (std::max_element(votes.begin(), votes.end()) - votes.begin()) - 1;
            // } else{
            //     best_direction = -1;
            // }
            // for (int m = 0; m < 4; m++) {
            //     std::cout << "direction " << m - 1 << " votes: " << votes[m] << std::endl;
            // }
            // std::cout << "best direction: " << best_direction << std::endl;
            // timer.end();
            // timer.print("Tuner");
        }

        int get_best_direction(){
            return best_direction;
        }

        void copy_in_level_0_info(T level_0_value_range, T level_0_error_bound_, std::vector<uint32_t> level_0_sizes_){
            value_range = level_0_value_range;
            level_0_error_bound = level_0_error_bound_;
            level_0_sizes = level_0_sizes_;
        }

        ~CoeffProfilingSamplingTuner(){}

        void print() const {
            std::cout << "Profiling Sampling Tuner with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Encoder:"; encoder.print();
            std::cout << "Interperter:"; interpreter.print();
        }
    private:
        void test_direct_BP_encoding(uint8_t num_bitplanes, uint32_t block_size=5){
            level_error_bounds.clear();
            // level_error_bounds.push_back(level_0_error_bound);
            level_sizes.clear();
            // level_sizes.push_back(level_0_sizes);
            size_t num_dims = dimensions.size();
            size_t block_buffer_size = 1;
            for(int i=0; i<num_dims; i++){
                block_buffer_size *= block_size;
            }
            size_t buffer_size = block_buffer_size * sampled_blocks.size();
            T* buffer = (T*) malloc(buffer_size * sizeof(T));
            T* buffer_ptr = buffer;
            for(int i=0; i<sampled_blocks.size(); i++){
                memcpy(buffer_ptr, sampled_blocks[i].data(), block_buffer_size * sizeof(T));
                buffer_ptr += block_buffer_size;
            }
            T level_max_error = compute_max_abs_value(buffer, buffer_size);
            if(negabinary) level_error_bounds.push_back(level_max_error * 4);
            else level_error_bounds.push_back(level_max_error);
            int level_exp = 0;
            frexp(level_max_error, &level_exp);
            std::vector<uint32_t> stream_sizes;
            auto streams = encoder.encode(buffer, buffer_size, level_exp, num_bitplanes, stream_sizes);
            uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
            level_sizes.push_back(stream_sizes);
        }
        void test_refactor(int direction, uint8_t target_level, uint8_t num_bitplanes, uint32_t block_size=5){
            level_error_bounds.clear();
            // level_error_bounds.push_back(level_0_error_bound);
            level_sizes.clear();
            // level_sizes.push_back(level_0_sizes);
            decomposer.direction = direction;
            size_t num_dims = dimensions.size();
            std::vector<uint32_t> block_dims;
            std::vector<uint32_t> block_dims_uint32;
            for(int i=0; i<num_dims; i++){
                block_dims.push_back(block_size);
                if(direction != i) block_dims_uint32.push_back(block_size);
            }
            auto level_dims = compute_level_dims_new(block_dims_uint32, target_level);
            auto level_elements = compute_level_elements(level_dims, target_level);
            for(int i=0; i<level_elements.size(); i++){
                level_elements[i] *= block_dims[direction];
            }
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
                auto block_level_buffers = decomposer.decompose_interleave_combine_levels(sampled_blocks_copied[b].data(), block_dims, target_level);
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

        uint32_t test_reconstruct(double tolerance, uint8_t target_level){
            std::vector<std::vector<double>> level_errors;
            MaxErrorCollector<T> collector = MaxErrorCollector<T>();
            for(int i=0; i<=target_level; i++){
                auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                level_errors.push_back(collected_error);
            }
            std::vector<uint8_t> level_num_bitplanes(target_level + 1, 0);
            // std::cout << "level_error_bounds: " << std::endl;
            // for(int i=0; i<= target_level; i++){
            //     std::cout << level_error_bounds[i] << " ";
            // }
            // std::cout << std::endl;
            auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance, level_num_bitplanes);
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

        Decomposer decomposer;
        Encoder encoder;
        Compressor compressor;
        SizeInterpreter interpreter;
        std::vector<std::vector<T>> sampled_blocks;
        std::vector<size_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<std::vector<uint32_t>> level_sizes;
        int best_direction = -1;
        T value_range;
        T level_0_error_bound;
        std::vector<uint32_t> level_0_sizes;
    public:
        bool negabinary = false;
    };
}
#endif