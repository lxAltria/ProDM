#ifndef _MDR_ORDERED_REFACTOR_NEW_HPP
#define _MDR_ORDERED_REFACTOR_NEW_HPP

#include "RefactorInterface.hpp"
#include "MDR/Decomposer/Decomposer.hpp"
#include "MDR/Interleaver/Interleaver.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/ErrorEstimator/ErrorEstimator.hpp"
#include "MDR/ErrorCollector/ErrorCollector.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/SizeInterpreter/SizeInterpreter.hpp"
#include "MDR/Writer/Writer.hpp"
#include "MDR/RefactorUtils.hpp"
#include <queue>

namespace MDR {

    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, error estimator and error collector
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class ErrorEstimator, class Writer>
    class OrderedRefactor_new : public concepts::RefactorInterface<T> {
    public:
        OrderedRefactor_new(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, ErrorEstimator error_estimator, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), collector(collector), error_estimator(error_estimator), writer(writer) {}

        // buffer: [metadata_size(uint32_t)][metadata][data]
        // return buffer size
        uint32_t refactor_to_buffer(T const * data_,
                                     const std::vector<uint32_t> &dims,
                                     uint8_t target_level,
                                     uint8_t num_bitplanes,
                                     uint8_t * buffer)
        {
            Timer timer;
            timer.start();
            dimensions = dims;
            uint32_t num_elements = 1;
            for (const auto &dim : dimensions)
            {
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);

            if (refactor(target_level, num_bitplanes))
            {
                timer.end();
                timer.print("Refactor");
            }

            std::vector<std::vector<double>> level_abs_errors;
            std::vector<std::vector<double>> &level_errors = level_squared_errors;
            if (std::is_base_of<MaxErrorEstimator<T>, ErrorEstimator>::value)
            {
                std::cout << "Computing absolute error" << std::endl;
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for (int i = 0; i <= target_level; i++)
                {
                    auto collected_error =
                        collector.collect_level_error(NULL, 0,
                                                      level_sizes[i].size(),
                                                      level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }
            else if (std::is_base_of<SquaredErrorEstimator<T>, ErrorEstimator>::value)
            {
                std::cout << "Using level squared error directly" << std::endl;
            }
            else
            {
                std::cerr << "Customized error estimator not supported yet"
                          << std::endl;
                exit(-1);
            }

            chunk_order = get_chunks_order(level_errors, error_perstep);

            // metadata
            uint32_t metadata_size = 0;
            uint8_t *metadata = get_metadata(metadata_size);

            // calculate data size
            std::vector<uint8_t> consumed(level_sizes.size(), 0);
            uint64_t total_data_size_64 = 0;
            for (uint8_t lev : chunk_order)
            {
                uint8_t j = consumed[lev];
                total_data_size_64 += level_sizes[lev][j];
                consumed[lev] = j + 1;
            }

            // 64 bit to 32 bit
            uint32_t data_size = static_cast<uint32_t>(total_data_size_64);

            // buffer size = sizeof(uint32_t) + metadata + data
            uint32_t buffer_size = static_cast<uint32_t>(sizeof(uint32_t)) +
                          metadata_size + data_size;

            // buffer shoule be already allocated
            uint8_t *p = buffer;
            if (!buffer) {
                std::cerr << "Buffer not allocated"
                          << std::endl;
                exit(-1);
            }

            // write metadata_size（uint32_t）
            std::memcpy(p, &metadata_size, sizeof(uint32_t));
            p += sizeof(uint32_t);

            // write metadata
            std::memcpy(p, metadata, metadata_size);
            p += metadata_size;

            // write compressed data by chunk_order
            std::fill(consumed.begin(), consumed.end(), 0);
            for (uint8_t lev : chunk_order)
            {
                uint8_t j = consumed[lev];
                uint32_t sz = level_sizes[lev][j];
                std::memcpy(p, level_components[lev][j], sz);
                p += sz;
                consumed[lev] = j + 1;
            }

            free(metadata);
            for (int i = 0; i < static_cast<int>(level_components.size()); i++)
            {
                for (int j = 0; j < static_cast<int>(level_components[i].size());
                     j++)
                {
                    free(level_components[i][j]);
                }
            }

            return buffer_size;
        }

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
            if(std::is_base_of<MaxErrorEstimator<T>, ErrorEstimator>::value){
                std::cout << "Computing absolute error" << std::endl;
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for(int i=0; i<level_sizes.size(); i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_sizes[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }
                level_errors = level_abs_errors;
            }
            else if(std::is_base_of<SquaredErrorEstimator<T>, ErrorEstimator>::value){
                std::cout << "Using level squared error directly" << std::endl;
            }
            else{
                std::cerr << "Customized error estimator not supported yet" << std::endl;
                exit(-1);
            }

            if(greedy) chunk_order = get_chunks_order(level_errors, error_perstep);
            if(bfs) {
                auto size_interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(error_estimator);
                std::vector<double> tmp_tolerances = {5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10, 5e-11, 1e-11, 0};
                std::vector<uint8_t> tmp_index(level_sizes.size(), 0);
                for(int j=0; j<tmp_tolerances.size(); j++){
                    tmp_tolerances[j] *= value_range;
                    auto output = size_interpreter.interpret_retrieve_size(level_sizes, level_errors, tmp_tolerances[j], tmp_index);
                }
                chunk_order = size_interpreter.get_path();
                error_perstep = size_interpreter.get_error_perstep();
            }
            write_metadata();

            // Calculate total size
            std::vector<uint8_t> consumed(level_sizes.size(), 0);
            uint64_t total_size_64 = 0;
            for (uint8_t lev : chunk_order) {
                uint8_t j = consumed[lev];
                total_size_64 += level_sizes[lev][j];
                consumed[lev] = j + 1;
            }

            // Attention: file size cannot be over 4GB
            uint32_t total_size = static_cast<uint32_t>(total_size_64);

            std::vector<uint8_t> packed(total_size);
            uint8_t* dst = packed.data();
            // copy every chunk
            std::fill(consumed.begin(), consumed.end(), 0);
            for (uint8_t lev : chunk_order) {
                uint8_t j = consumed[lev];
                uint32_t sz = level_sizes[lev][j];
                std::memcpy(dst, level_components[lev][j], sz);
                dst += sz;
                consumed[lev] = j + 1;
            }

            writer.write_components(packed.data(), total_size);
            
            // level_num.clear();
            // level_num.push_back(1);
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
                + get_size(level_sizes)     
                + get_size(level_elements)                       
                + get_size(stopping_indices)
                + sizeof(uint8_t)    // negabinary
                + sizeof(uint16_t)   // chunk_num
                + get_size(chunk_order)                                     
                + get_size(error_perstep);

            uint8_t* metadata = static_cast<uint8_t*>(malloc(metadata_size));
            uint8_t* p = metadata;

            *(p++) = dimensions.size();
            serialize(dimensions, p);
            *(p++) = level_error_bounds.size();
            serialize(level_error_bounds, p);
            serialize(level_sizes, p);
            serialize(level_elements, p);
            serialize(stopping_indices, p);
            *(p++) = static_cast<uint8_t>(negabinary);

            const uint16_t chunk_num = chunk_order.size();
            memcpy(p, &chunk_num, sizeof(uint16_t));
            p += sizeof(uint16_t);
            serialize(chunk_order, p);
            serialize(error_perstep, p);

            return metadata;
        }

        void write_metadata() const {
            uint32_t metadata_size;
            uint8_t* metadata = get_metadata(metadata_size);

            writer.write_metadata(metadata, metadata_size);
            free(metadata);
        }

        ~OrderedRefactor_new(){}

        void print() const {
            std::cout << "Composed refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        std::vector<uint8_t> get_chunks_order(const std::vector<std::vector<double>>& level_errors, std::vector<double>& error_perstep) const {
            // for(int i=0; i<level_errors.size(); i++){
            //     for(int j=0; j<level_errors[i].size(); j++){
            //         std::cout << level_errors[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            std::vector<uint8_t> index(level_sizes.size(), 0);
            int num_levels = level_sizes.size();
            std::vector<uint32_t> retrieve_sizes(num_levels, 0);
            double accumulated_error = 0;
            for(int i=0; i<num_levels; i++){
                accumulated_error += error_estimator.estimate_error(level_errors[i][index[i]], i, num_levels);
            }
            std::priority_queue<UnitErrorGain, std::vector<UnitErrorGain>, CompareUnitErrorGain> heap;
            // identify minimal level
            double min_error = accumulated_error;

            std::vector<uint8_t> order;
            order.reserve([&](){
                size_t s=0; for(const auto& v: level_sizes) s+=v.size(); return s;
            }()); // allocate enough space for it

            for(int i=0; i<num_levels; i++){
                min_error -= error_estimator.estimate_error(level_errors[i][index[i]], i, num_levels);
                min_error += error_estimator.estimate_error(level_errors[i].back(), i, num_levels);
                // fetch the first component if index is 0
                if(index[i] == 0){
                    retrieve_sizes[i] += level_sizes[i][index[i]];
                    accumulated_error -= error_estimator.estimate_error(level_errors[i][index[i]], i, num_levels);
                    accumulated_error += error_estimator.estimate_error(level_errors[i][index[i] + 1], i, num_levels);
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
                retrieve_sizes[i] += level_sizes[i][j];
                accumulated_error -= error_estimator.estimate_error(level_errors[i][j], i, num_levels);
                accumulated_error += error_estimator.estimate_error(level_errors[i][j + 1], i, num_levels);
                index[i] ++;
                if(index[i] != level_sizes[i].size()){
                    double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i] + 1], i);
                    heap.push(UnitErrorGain(error_gain / level_sizes[i][index[i]], i));
                }
                order.push_back(static_cast<uint8_t>(i));
                error_perstep.push_back(accumulated_error);
            }
            return order;
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
            std::vector<std::vector<T>> level_buffers = decomposer.decompose_interleave(data.data(), dimensions, target_level);
            // MGARD::writefile("decomposed_coeff.dat", data.data(), data.size());
            // timer.end();
            // timer.print("Decompose");

            // encode level by level
            level_error_bounds.clear();
            level_squared_errors.clear();
            level_components.clear();
            level_sizes.clear();
            // auto level_dims = compute_level_dims(dimensions, target_level);
            // auto level_elements = compute_level_elements(level_dims, target_level);
            // std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
            SquaredErrorCollector<T> s_collector = SquaredErrorCollector<T>();
            for(int i=0; i<level_buffers.size(); i++){
                // timer.start();
                // const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                // T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
                // extract level i component
                // interleaver.interleave(data.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
                // compute max coefficient as level error bound
                level_elements.push_back(level_buffers[i].size());
                T level_max_error = compute_max_abs_value(level_buffers[i].data(), level_elements[i]);
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
                auto streams = encoder.encode(level_buffers[i].data(), level_elements[i], level_exp, num_bitplanes, stream_sizes);
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
        ErrorEstimator error_estimator;
        Writer writer;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<uint8_t> chunk_order;
        std::vector<double> error_perstep;
        std::vector<uint32_t> level_elements;
    public:
        bool negabinary = false;
        bool greedy = false;
        bool bfs = false;
    };
}
#endif