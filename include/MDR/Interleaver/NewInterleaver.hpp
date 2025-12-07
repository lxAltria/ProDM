#ifndef _MDR_DIRECT_INTERLEAVER_NEW_HPP
#define _MDR_DIRECT_INTERLEAVER_NEW_HPP

#include "InterleaverInterface.hpp"
#include <cassert>

namespace MDR {
    // direct interleaver with in-order recording
    template<class T>
    class DirectInterleaver_new : public concepts::InterleaverInterface<T> {
    public:
        DirectInterleaver_new(){}
        void interleave(T const * data, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * buffer, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {

        }
        void interleave(T const * data, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * buffer, int level, int target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            int h = 1 << (target_level - level);
            // std::cout << "target_level: " << target_level << ", level: " << level << ", h: " << h << std::endl;
            if(dims.size() == 1){
                if(level){
                    size_t begin = 0;
                    size_t end = dims[0] - 1;
                    uint32_t count = interleave_1D(begin, end, h, buffer, const_cast<T*>(data));
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0] - dims_coasre[0]));
                }
                else{
                    size_t begin = 0;
                    size_t end = dims[0] - 1;
                    uint32_t count = interleave_1D_level_0(begin, end, h, buffer, const_cast<T*>(data));
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0] - dims_coasre[0]));
                }
            }
            else if(dims.size() == 2){
                if(level){
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    uint32_t count = interleave_2D(buffer, const_cast<T*>(data), n1, n2, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1] - dims_coasre[0]*dims_coasre[1]));
                }
                else{
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    uint32_t count = interleave_2D_level_0(buffer, const_cast<T*>(data), n1, n2, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1] - dims_coasre[0]*dims_coasre[1]));
                }
            }
            else if(dims.size() == 3){  
                if(level){
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    size_t n3 = dims[2];
                    uint32_t count = interleave_3D(buffer, const_cast<T*>(data), n1, n2, n3, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1]*dims_fine[2] - dims_coasre[0]*dims_coasre[1]*dims_coasre[2]));
                }
                else{
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    size_t n3 = dims[2];
                    uint32_t count = interleave_3D_level_0(buffer, const_cast<T*>(data), n1, n2, n3, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1]*dims_fine[2] - dims_coasre[0]*dims_coasre[1]*dims_coasre[2]));
                }
            }
            else{
                std::cout << "Dimension higher than 4 is not supported\n";
                exit(-1);
            }
        }
        void reposition(T const * buffer, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * data, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {

        }
        void reposition(T const * buffer, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * data, int level, int target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            int h = 1 << (target_level - level);
            // std::cout << "target_level: " << target_level << ", level: " << level << ", h: " << h << std::endl;
            if(dims.size() == 1){
                if(level){
                    size_t begin = 0;
                    size_t end = dims[0] - 1;
                    uint32_t count = reposition_1D(begin, end, h, const_cast<T*>(buffer), data);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0] - dims_coasre[0]));
                }
                else{
                    size_t begin = 0;
                    size_t end = dims[0] - 1;
                    uint32_t count = reposition_1D_level_0(begin, end, h, const_cast<T*>(buffer), data);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0] - dims_coasre[0]));
                }
            }
            else if(dims.size() == 2){
                if(level){
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    uint32_t count = reposition_2D(const_cast<T*>(buffer), data, n1, n2, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1] - dims_coasre[0]*dims_coasre[1]));
                }
                else{
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    uint32_t count = reposition_2D_level_0(const_cast<T*>(buffer), data, n1, n2, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1] - dims_coasre[0]*dims_coasre[1]));
                }
            }
            else if(dims.size() == 3){
                if(level){
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    size_t n3 = dims[2];
                    uint32_t count = reposition_3D(const_cast<T*>(buffer), data, n1, n2, n3, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1]*dims_fine[2] - dims_coasre[0]*dims_coasre[1]*dims_coasre[2]));
                }
                else{
                    size_t n1 = dims[0];
                    size_t n2 = dims[1];
                    size_t n3 = dims[2];
                    uint32_t count = reposition_3D_level_0(const_cast<T*>(buffer), data, n1, n2, n3, h);
                    // std::cout << "count: " << count << std::endl;
                    assert(count == (dims_fine[0]*dims_fine[1]*dims_fine[2] - dims_coasre[0]*dims_coasre[1]*dims_coasre[2]));
                }
            }
            else{
                std::cout << "Dimension higher than 4 is not supported\n";
                exit(-1);
            }
        }
        void print() const {
            std::cout << "Direct interleaver new" << std::endl;
        }
    private:
        uint32_t interleave_1D_level_0(const size_t begin, const size_t end, const size_t stride, T * buffer, T * data){
            uint32_t count = 0;
            for(size_t i=begin; i<=end; i+=stride){
                buffer[count++] = data[i];
            }
            return count;
        }
        uint32_t interleave_2D_level_0(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t h){
            uint32_t count_2D = 0;
            uint32_t count;
            size_t stride_n2 = h;
            size_t stride_n1 = n2 * h;

            T * cur_data_pos = data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n2_begin = 0;
            size_t n2_end = n2 - 1;
            for(size_t i=0; i<n1; i+=h){
                count = interleave_1D_level_0(n2_begin, n2_end, stride_n2, cur_buffer_pos, cur_data_pos);
                cur_data_pos += stride_n1;
                cur_buffer_pos += count;
                count_2D += count;
            }
            return count_2D;
        }
        uint32_t interleave_3D_level_0(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
            uint32_t count_3D = 0;
            uint32_t count;
            size_t stride_n3 = h;
            size_t stride_n2 = n3*h;
            size_t stride_n1 = n2 * n3 * h;

            T * cur_data_pos = data_pos;
            T * temp_data_pos = cur_data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n3_begin = 0;
            size_t n3_end = n3 - 1;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n2; j+=h){
                    count = interleave_1D_level_0(n3_begin, n3_end, h, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += n3*h;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            return count_3D;
        }
        uint32_t interleave_1D(const size_t begin, const size_t end, const size_t stride, T * buffer, T * data){
            uint32_t count = 0;
            size_t n = (end - begin) / stride + 1;
            size_t i = 1;
            for(i=1; i+1<n; i+=2){
                size_t c = begin + i * stride;
                buffer[count ++] = data[c];
            }
            if(n % 2 == 0){
                size_t c = begin + (n - 1) * stride;
                buffer[count ++] = data[c];
            }
            return count;
        }
        uint32_t interleave_2D(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t h){
            uint32_t count2D = 0;
            uint32_t count;
            size_t stride_n1 = n2 * h;
            size_t stride_n2 = h;

            T * cur_data_pos = data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n2_begin = 0;
            size_t n2_end = n2 - 1;
            for(size_t i=0; i<n1; i+=h){
                count = interleave_1D(n2_begin, n2_end, stride_n2, cur_buffer_pos, cur_data_pos);
                cur_data_pos += n2 * h;
                cur_buffer_pos += count;
                count2D += count;
            }
            // compute vertical difference
            cur_data_pos = data_pos;
            h = h << 1;
            size_t n1_begin = 0;
            size_t n1_end = (n1 - 1) * n2;
            for(size_t i=0; i<n2; i+=h){
                count = interleave_1D(n1_begin, n1_end, stride_n1, cur_buffer_pos, cur_data_pos);
                cur_data_pos += h;
                cur_buffer_pos += count;
                count2D += count;
            }
            return count2D;
        }
        uint32_t interleave_3D(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
            uint32_t count_3D = 0;
            uint32_t count;
            size_t stride_n1 = n2 * n3 * h;
            size_t stride_n2 = n3 * h;

            T * cur_data_pos = data_pos;
            T * temp_data_pos = cur_data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n3_begin = 0;
            size_t n3_end = n3 - 1;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n2; j+=h){
                    count = interleave_1D(n3_begin, n3_end, h, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += n3*h;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            size_t h2x = h << 1;
            cur_data_pos = data_pos;
            size_t n2_begin = 0;
            size_t n2_end = (n2-1) * n3;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n3; j+=h2x){
                    count = interleave_1D(n2_begin, n2_end, stride_n2, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += h2x;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            cur_data_pos = data_pos;
            size_t n1_begin = 0;
            size_t n1_end = (n1-1) * n2 * n3;
            for(size_t i=0; i<n2; i+=h2x){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n3; j+=h2x){
                    count = interleave_1D(n1_begin, n1_end, stride_n1, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += h2x;
                    count_3D += count;
                }
                cur_data_pos += n3 * h2x;
            }

            return count_3D;
        }
        uint32_t reposition_1D_level_0(const size_t begin, const size_t end, const size_t stride, T * buffer, T * data){
            uint32_t count = 0;
            for(size_t i=begin; i<=end; i+=stride){
                data[i] = buffer[count++];
            }
            return count;
        }
        uint32_t reposition_2D_level_0(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t h){
            uint32_t count_2D = 0;
            uint32_t count;
            size_t stride_n2 = h;
            size_t stride_n1 = n2 * h;

            T * cur_data_pos = data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n2_begin = 0;
            size_t n2_end = n2 - 1;
            for(size_t i=0; i<n1; i+=h){
                count = reposition_1D_level_0(n2_begin, n2_end, stride_n2, cur_buffer_pos, cur_data_pos);
                cur_data_pos += stride_n1;
                cur_buffer_pos += count;
                count_2D += count;
            }
            return count_2D;
        }
        uint32_t reposition_3D_level_0(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
            uint32_t count_3D = 0;
            uint32_t count;
            size_t stride_n3 = h;
            size_t stride_n1 = n2 * n3 * h;

            T * cur_data_pos = data_pos;
            T * temp_data_pos = cur_data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n3_begin = 0;
            size_t n3_end = n3 - 1;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n2; j+=h){
                    count = reposition_1D_level_0(n3_begin, n3_end, h, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += n3*h;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            return count_3D;
        }
        uint32_t reposition_1D(const size_t begin, const size_t end, const size_t stride, T * buffer, T * data){
            uint32_t count = 0;
            size_t n = (end - begin) / stride + 1;
            size_t i = 1;
            for(i=1; i+1<n; i+=2){
                size_t c = begin + i * stride;
                data[c] = buffer[count++];
            }
            if(n % 2 == 0){
                size_t c = begin + (n - 1) * stride;
                data[c] = buffer[count++];
            }
            return count;
        }
        uint32_t reposition_2D(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t h){
            uint32_t count2D = 0;
            uint32_t count;
            size_t stride_n1 = n2 * h;
            size_t stride_n2 = h;

            T * cur_data_pos = data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n2_begin = 0;
            size_t n2_end = n2 - 1;
            for(size_t i=0; i<n1; i+=h){
                count = reposition_1D(n2_begin, n2_end, stride_n2, cur_buffer_pos, cur_data_pos);
                cur_data_pos += n2 * h;
                cur_buffer_pos += count;
                count2D += count;
            }
            // compute vertical difference
            cur_data_pos = data_pos;
            h = h << 1;
            size_t n1_begin = 0;
            size_t n1_end = (n1 - 1) * n2;
            for(size_t i=0; i<n2; i+=h){
                count = reposition_1D(n1_begin, n1_end, stride_n1, cur_buffer_pos, cur_data_pos);
                cur_data_pos += h;
                cur_buffer_pos += count;
                count2D += count;
            }
            return count2D;
        }
        uint32_t reposition_3D(T * buffer_pos, T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
            uint32_t count_3D = 0;
            uint32_t count;
            size_t stride_n1 = n2 * n3 * h;
            size_t stride_n2 = n3 * h;

            T * cur_data_pos = data_pos;
            T * temp_data_pos = cur_data_pos;
            T * cur_buffer_pos = buffer_pos;
            size_t n3_begin = 0;
            size_t n3_end = n3 - 1;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n2; j+=h){
                    count = reposition_1D(n3_begin, n3_end, h, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += n3*h;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            size_t h2x = h << 1;
            cur_data_pos = data_pos;
            size_t n2_begin = 0;
            size_t n2_end = (n2-1) * n3;
            for(size_t i=0; i<n1; i+=h){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n3; j+=h2x){
                    count = reposition_1D(n2_begin, n2_end, stride_n2, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += h2x;
                    count_3D += count;
                }
                cur_data_pos += n2*n3*h;
            }

            cur_data_pos = data_pos;
            size_t n1_begin = 0;
            size_t n1_end = (n1-1) * n2 * n3;
            for(size_t i=0; i<n2; i+=h2x){
                temp_data_pos = cur_data_pos;
                for(size_t j=0; j<n3; j+=h2x){
                    count = reposition_1D(n1_begin, n1_end, stride_n1, cur_buffer_pos, temp_data_pos);
                    cur_buffer_pos += count;
                    temp_data_pos += h2x;
                    count_3D += count;
                }
                cur_data_pos += n3 * h2x;
            }
            return count_3D;
        }
    };
}
#endif
