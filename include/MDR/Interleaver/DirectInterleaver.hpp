#ifndef _MDR_DIRECT_INTERLEAVER_HPP
#define _MDR_DIRECT_INTERLEAVER_HPP

#include "InterleaverInterface.hpp"

namespace MDR {
    // direct interleaver with in-order recording
    template<class T>
    class DirectInterleaver : public concepts::InterleaverInterface<T> {
    public:
        DirectInterleaver(){}
        void interleave(T const * data, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * buffer, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            if(dims.size() == 1){
                uint32_t count = 0;
                for(int i=dims_coasre[0]; i<dims_fine[0]; i++){
                    buffer[count ++] = data[i];
                }
            }
            else if(dims.size() == 3){
                uint32_t dim0_offset = strides.size() ? strides[0] : (dims[1] * dims[2]);
                uint32_t dim1_offset = strides.size() ? strides[1] : dims[2];
                uint32_t count = 0;
                for(int i=0; i<dims_fine[0]; i++){
                    for(int j=0; j<dims_fine[1]; j++){
                        for(int k=0; k<dims_fine[2]; k++){
                            if((i < dims_coasre[0]) && (j < dims_coasre[1]) && (k < dims_coasre[2]))
                                continue;
                            buffer[count ++] = data[i*dim0_offset + j*dim1_offset + k];
                        }
                    }
                }                
            }
            else{
                std::cout << "Dimension higher than 4 is not supported\n";
                exit(-1);
            }
        }
        void reposition(T const * buffer, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * data, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            if(dims.size() == 1){
                uint32_t count = 0;
                for(int i=dims_coasre[0]; i<dims_fine[0]; i++){
                    data[i] = buffer[count ++];
                }
            }
            else if(dims.size() == 3){
                uint32_t dim0_offset = strides.size() ? strides[0] : (dims[1] * dims[2]);
                uint32_t dim1_offset = strides.size() ? strides[1] : dims[2];
                uint32_t count = 0;
                for(int i=0; i<dims_fine[0]; i++){
                    for(int j=0; j<dims_fine[1]; j++){
                        for(int k=0; k<dims_fine[2]; k++){
                            if((i < dims_coasre[0]) && (j < dims_coasre[1]) && (k < dims_coasre[2]))
                                continue;
                            data[i*dim0_offset + j*dim1_offset + k] = buffer[count ++];
                        }
                    }
                }
            }
            else{
                std::cout << "Dimension higher than 4 is not supported\n";
                exit(-1);
            }
        }
        void print() const {
            std::cout << "Direct interleaver" << std::endl;
        }
    };
}
#endif
