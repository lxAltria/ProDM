#ifndef _PDR_APPROXIMATION_BASED_DELTA_REFACTOR_HPP
#define _PDR_APPROXIMATION_BASED_DELTA_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "PDR/Approximator/Approximator.hpp"

namespace PDR {

    // an approximation-based scientific data refactor: compose an approximation algorithm, encoder, and lossless compressor
    template<class T, class Approximator>
    class ApproximationBasedDeltaRefactor : public concepts::RefactorInterface<T> {
    public:
        ApproximationBasedDeltaRefactor(Approximator approximator)
            : approximator(approximator) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){

            num_segments = std::is_same<T, float>::value ? 7 : 16;
            Timer timer;
            timer.start();
            dimensions = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            double approximator_eb = 0.1;
            T max_val = data[0];
            T min_val = data[0];
            for(int i=1; i<num_elements; i++){
                if(data[i] > max_val) max_val = data[i];
                if(data[i] < min_val) min_val = data[i];
            }
            value_range = max_val - min_val;
            approximator_eb *= value_range;
            // refactor
            level_sizes.clear();
            for(int i=0; i<num_segments; i++){
                std::string filename = file_prefix + std::to_string(i);
                auto size = approximator.refactor_approximate(data.data(), dimensions, approximator_eb, filename);
                level_sizes.push_back(size);
                approximator_eb *= 0.1;
            }
            timer.end();
            timer.print("Refactor");

            write_metadata();
        }

        void write_metadata() const {
            uint32_t metadata_size = sizeof(uint8_t) + get_size(dimensions) // dimensions
                            + sizeof(double) + sizeof(int) + get_size(level_sizes); // one byte for whether negabinary encoding is used 
            uint8_t * metadata = (uint8_t *) malloc(metadata_size);
            uint8_t * metadata_pos = metadata;
            *(metadata_pos ++) = (uint8_t) dimensions.size();
            serialize(dimensions, metadata_pos);
            *reinterpret_cast<double*>(metadata_pos) = value_range;
            metadata_pos += sizeof(double);
            *reinterpret_cast<int*>(metadata_pos) = num_segments;            
            metadata_pos += sizeof(int);
            serialize(level_sizes, metadata_pos);
            FILE * file = fopen(metadata_file.c_str(), "w");
            fwrite(metadata, 1, metadata_size, file);
            fclose(file);
            free(metadata);
        }

        ~ApproximationBasedDeltaRefactor(){}

        void print() const {
            std::cout << "Approximation-based delta refactor with the following components." << std::endl;
            std::cout << "Approximator: "; approximator.print();
        }

    private:

        Approximator approximator;
        int num_segments = 0;
        T value_range = 0;
        std::string metadata_file = "refactored_data/metadata.bin";
        std::string file_prefix = "refactored_data/delta_segment_";
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<uint32_t> level_sizes;
    };
}
#endif

