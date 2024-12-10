#ifndef _PDR_COMPOSED_DELTA_RECONSTRUCTOR_HPP
#define _PDR_COMPOSED_DELTA_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "PDR/Approximator/Approximator.hpp"

namespace PDR {
    // an approximation-based scientific data reconstructor: inverse operator of approximation-based refactor
    template<class T, class Approximator>
    class ApproximationBasedDeltaReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        ApproximationBasedDeltaReconstructor(Approximator approximator)
            : approximator(approximator){}

        T * reconstruct(double tolerance){
            return progressive_reconstruct(tolerance);
        }

        T * progressive_reconstruct(double tolerance){
            std::vector<T> buffer(num_elements);
            while(approximator_eb > tolerance){
                if(current_segment >= num_segments) break;
                // fetch one more segment
                std::string filename = file_prefix + std::to_string(current_segment);
                approximator.reconstruct_approximate(buffer.data(), dimensions, filename);
                for(int i=0; i<num_elements; i++){
                    data[i] += buffer[i];
                }
                retrieval_size += level_sizes[current_segment];
                current_segment ++;
                approximator_eb *= 0.1;
            }
            return data.data();
        }

        T * progressive_reconstruct(double tolerance, int max_level){
            return progressive_reconstruct(tolerance);
        }

        void load_metadata(){
            FILE * file = fopen(metadata_file.c_str(), "r");
            fseek(file, 0, SEEK_END);
            uint32_t num_bytes = ftell(file);
            rewind(file);
            uint8_t * metadata = (uint8_t *) malloc(num_bytes);
            fread(metadata, 1, num_bytes, file);
            fclose(file);            
            uint8_t const * metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos ++);
            deserialize(metadata_pos, num_dims, dimensions);
            approximator_eb = *reinterpret_cast<const double*>(metadata_pos);
            metadata_pos += sizeof(double);
            num_segments = *reinterpret_cast<const int*>(metadata_pos);
            metadata_pos += sizeof(int);
            deserialize(metadata_pos, num_segments, level_sizes);
            strides = std::vector<uint32_t>(dimensions.size());
            uint32_t stride = 1;
            for(int i=dimensions.size()-1; i>=0; i--){
                strides[i] = stride;
                stride *= dimensions[i];
            }
            data = std::vector<T>(stride, 0);
            num_elements = stride;
            free(metadata);
        }

        const std::vector<uint32_t>& get_dimensions(){
            return dimensions;
        }

        size_t get_retrieved_size(){
            return retrieval_size;
        }

        ~ApproximationBasedDeltaReconstructor(){}

        void print() const {
            std::cout << "Approximation-based delta reconstructor with the following components." << std::endl;
            std::cout << "Approximatorr: "; approximator.print();
        }

    private:
        Approximator approximator;
        double approximator_eb = 0;
        int num_segments = 0;
        int current_segment = 0;
        size_t retrieval_size = 0;
        size_t num_elements = 0;
        std::string metadata_file = "refactored_data/metadata.bin";
        std::string file_prefix = "refactored_data/delta_segment_";
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<uint32_t> level_sizes;
        std::vector<uint32_t> strides;
    };
}
#endif

