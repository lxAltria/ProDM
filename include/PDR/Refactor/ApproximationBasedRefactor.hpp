#ifndef _PDR_APPROXIMATION_BASED_REFACTOR_HPP
#define _PDR_APPROXIMATION_BASED_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "PDR/Approximator/Approximator.hpp"
#include "MDR/BitplaneEncoder/BitplaneEncoder.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"
#include "MDR/Writer/Writer.hpp"
#include "MDR/RefactorUtils.hpp"

using namespace MDR;

namespace PDR {

    // an approximation-based scientific data refactor: compose an approximation algorithm, encoder, and lossless compressor
    template<class T, class Approximator, class Encoder, class Compressor, class Writer>
    class ApproximationBasedRefactor : public concepts::RefactorInterface<T> {
    public:
        ApproximationBasedRefactor(Approximator approximator, Encoder encoder, Compressor compressor, Writer writer)
            : approximator(approximator), encoder(encoder), compressor(compressor), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){
            // Timer timer;
            // timer.start();
            dimensions = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            // if refactor successfully
            if(refactor(num_bitplanes)){
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
                            + sizeof(size_t)
                            + sizeof(uint8_t) + get_size(level_error_bounds) 
                            + get_size(level_sizes) // level information
                            + get_size(stopping_indices) + get_size(level_num) + 1 + sizeof(T); // one byte for whether negabinary encoding is used 
            uint8_t * metadata = (uint8_t *) malloc(metadata_size);
            uint8_t * metadata_pos = metadata;
            *(metadata_pos ++) = (uint8_t) dimensions.size();
            serialize(dimensions, metadata_pos);
            *reinterpret_cast<size_t*>(metadata_pos) = approximator_size;
            metadata_pos += sizeof(size_t);
            *(metadata_pos ++) = (uint8_t) 1; // level = 1
            serialize(level_error_bounds, metadata_pos);
            serialize(level_sizes, metadata_pos);
            serialize(stopping_indices, metadata_pos);
            serialize(level_num, metadata_pos);
            *(metadata_pos ++) = (uint8_t) negabinary;
            *reinterpret_cast<T*>(metadata_pos) = approximator_eb;
            writer.write_metadata(metadata, metadata_size);
            free(metadata);
        }

        ~ApproximationBasedRefactor(){}

        void print() const {
            std::cout << "Approximation-based refactor with the following components." << std::endl;
            std::cout << "Approximator: "; approximator.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        bool refactor(uint8_t num_bitplanes){

            auto num_elements = data.size();
            T max_val = data[0];
            T min_val = data[0];
            for(int i=1; i<num_elements; i++){
                if(data[i] > max_val) max_val = data[i];
                if(data[i] < min_val) min_val = data[i];
            }
            approximator_eb *= (max_val - min_val);
            std::string approximator_path = writer.get_directory() + "approximator.dat";
            approximator_size = approximator.refactor_approximate(data.data(), dimensions, approximator_eb, approximator_path);

            level_error_bounds.clear();
            level_components.clear();
            level_sizes.clear();
            
            if (mask.empty()){
                T level_max_error = compute_max_abs_value(data.data(), num_elements);
                // std::cout << "level_max_error = " << level_max_error << std::endl;

                // encoding
                if(negabinary) level_error_bounds.push_back(level_max_error * 4);
                else level_error_bounds.push_back(level_max_error);
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(data.data(), num_elements, level_exp, num_bitplanes, stream_sizes);

                // lossless
                uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                stopping_indices.push_back(stopping_index);
                level_components.push_back(streams);
                level_sizes.push_back(stream_sizes);
            }
            else{
                int num_valid_data = std::accumulate(mask.begin(), mask.end(), 0);
                std::vector<T> filtered_data(num_valid_data);
                int filtered_index = 0;
                for(int i=0; i<mask.size(); i++){
                    if(mask[i]){
                        filtered_data[filtered_index++]=data[i];
                    }
                }
                T level_max_error = compute_max_abs_value(filtered_data.data(), num_valid_data);
                // std::cout << "level_max_error = " << level_max_error << std::endl;

                // encoding
                if(negabinary) level_error_bounds.push_back(level_max_error * 4);
                else level_error_bounds.push_back(level_max_error);
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(filtered_data.data(), num_valid_data, level_exp, num_bitplanes, stream_sizes);

                // lossless
                uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                stopping_indices.push_back(stopping_index);
                level_components.push_back(streams);
                level_sizes.push_back(stream_sizes);
            }

            return true;
        }

        Approximator approximator;
        Encoder encoder;
        Compressor compressor;
        Writer writer;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        T approximator_eb = 0.001;
        size_t approximator_size = 0;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<uint32_t>> level_sizes;
    public:
        bool negabinary = false;
        std::vector<unsigned char> mask;
    };
}
#endif

