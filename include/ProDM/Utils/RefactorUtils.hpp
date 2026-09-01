#ifndef _MDR_REFACTOR_UTILS_HPP
#define _MDR_REFACTOR_UTILS_HPP

#include <cassert>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <iostream>

namespace MDR {

    // MDR utility functions

    // MGARD related
    // TODO: put API in MGARD

    // compute level dimensions
    /*
        @params dims: input dimensions
        @params target_level: the target decomposition level
    */
    std::vector<std::vector<uint32_t>> compute_level_dims(const std::vector<uint32_t>& dims, uint32_t target_level){
        std::vector<std::vector<uint32_t>> level_dims;
        for(int i=0; i<=target_level; i++){
            level_dims.push_back(std::vector<uint32_t>(dims.size()));
        }
        for(int i=0; i<dims.size(); i++){
            int n = dims[i];
            for(int j=0; j<=target_level; j++){
                level_dims[target_level - j][i] = n;
                n = (n >> 1) + 1;
            }
        }
        return level_dims;
    }
    std::vector<std::vector<uint32_t>> compute_level_dims_new(const std::vector<uint32_t>& dims, uint32_t target_level){
        std::vector<std::vector<uint32_t>> level_dims;
        for(int i=0; i<=target_level; i++){
            level_dims.push_back(std::vector<uint32_t>(dims.size()));
        }
        for(int i=0; i<dims.size(); i++){
            int n = dims[i];
            for(int j=0; j<=target_level; j++){
                level_dims[target_level - j][i] = n;
                n = (n >> 1) + (n & 1);
            }
        }
        return level_dims;
    }

    // compute level elements
    /*
        @params level_dims: dimensions for all levels
        @params target_level: the target decomposition level
    */
    std::vector<uint32_t> compute_level_elements(const std::vector<std::vector<uint32_t>>& level_dims, int target_level){
        assert(level_dims.size());
        uint8_t num_dims = level_dims[0].size();
        std::vector<uint32_t> level_elements(level_dims.size());
        level_elements[0] = 1;
        for(int j=0; j<num_dims; j++){
            level_elements[0] *= level_dims[0][j];
        }
        uint32_t pre_num_elements = level_elements[0];
        for(int i=1; i<=target_level; i++){
            uint32_t num_elements = 1;
            for(int j=0; j<num_dims; j++){
                num_elements *= level_dims[i][j];
            }
            level_elements[i] = num_elements - pre_num_elements;
            pre_num_elements = num_elements;
        }
        return level_elements;
    }

    std::vector<uint32_t> compute_level_buffers_size(const std::vector<std::vector<uint32_t>>& level_dims, int target_level){
        assert(level_dims.size());
        size_t num_dims = level_dims[0].size();
        // size_t count;
        size_t size;
        std::vector<uint32_t> level_sizes(1 + num_dims * target_level);
        level_sizes[0] = 1;
        // std::cout << "level 0, dims = ";
        for(int i=0; i<num_dims; i++){
            level_sizes[0] *= level_dims[0][i];
            // std::cout << level_dims[0][i] << " ";
        }
        // std::cout << std::endl;
        for(size_t l=1; l<=target_level; l++){
            // count = 0;
            switch (num_dims){
                case 1:
                {
                    size = level_dims[l][0];
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    // assert(count == (level_dims[l][0] - level_dims[l-1][0]));
                    break;
                }
                case 2:
                {
                    size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1];
                    // std::cout << "level " << (l-1) * num_dims + 1 << ", dims = " << level_dims[l][0] - level_dims[l-1][0] << " " << level_dims[l-1][1] << std::endl;
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]);
                    // std::cout << "level " << (l-1) * num_dims + 2 << ", dims = " << level_dims[l][0] << " " << level_dims[l][1] - level_dims[l-1][1] << std::endl;
                    level_sizes[(l-1) * num_dims + 2] = size;
                    // count += size;
                    // assert(count == (level_dims[l][0]*level_dims[l][1] - level_dims[l-1][0]*level_dims[l-1][1]));
                    break;
                }
                case 3:
                {
                    // n1 (cur_n1 - pre_n1) * pre_n2 * pre_n3
                    size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1] * level_dims[l-1][2];
                    // std::cout << "level " << (l-1) * num_dims + 1 << ", dims = " << (level_dims[l][0] - level_dims[l-1][0]) << " " << level_dims[l-1][1] << " " << level_dims[l-1][2] << std::endl;
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    // n2 cur_n1 * (cur_n2 - pre_n2) * pre_n3
                    size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]) * level_dims[l-1][2];
                    // std::cout << "level " << (l-1) * num_dims + 2 << ", dims = " << level_dims[l][0] << " " << (level_dims[l][1] - level_dims[l-1][1]) << " " << level_dims[l-1][2] << std::endl;
                    level_sizes[(l-1) * num_dims + 2] = size;
                    // count += size;
                    // n3 cur_n1 * cur_n2 * (cur_n3 - pre_n3)
                    size = level_dims[l][0] * level_dims[l][1] * (level_dims[l][2] - level_dims[l-1][2]);
                    // std::cout << "level " << (l-1) * num_dims + 3 << ", dims = " << level_dims[l][0] << " " << level_dims[l][1] << " " << (level_dims[l][2] - level_dims[l-1][2]) << std::endl;
                    level_sizes[(l-1) * num_dims + 3] = size;
                    // count += size;
                    // assert(count == (level_dims[l][0]*level_dims[l][1]*level_dims[l][2] - level_dims[l-1][0]*level_dims[l-1][1]*level_dims[l-1][2]));
                    break;
                }
                default:
                    std::cerr << num_dims << "-Dimentional decomposition not implemented." << std::endl;
                    exit(-1);
            }
        }
        return level_sizes;
    }

    // Simple utility functions

    // compute maximum value in level
    /*
    @params data: level data
    @params n: number of level data points
    */
    template <class T>
    T compute_max_abs_value(const T * data, uint32_t n){
        T max_val = 0;
        for(int i=0; i<n; i++){
            T val = fabs(data[i]);
            if(val > max_val) max_val = val;
        }
        return max_val;
    }

    // Get size of vector
    template <class T>
    inline uint32_t get_size(const std::vector<T>& vec){
        return vec.size() * sizeof(T);
    }
    template <class T>
    uint32_t get_size(const std::vector<std::vector<T>>& vec){
        uint32_t size = 0;
        for(int i=0; i<vec.size(); i++){
            size += sizeof(uint32_t) + vec[i].size() * sizeof(T);
        }
        return size;
    }

    // Serialize/deserialize vectors
    // Auto-increment buffer position
    template <class T>
    inline void serialize(const std::vector<T>& vec, uint8_t *& buffer_pos){
        memcpy(buffer_pos, vec.data(), vec.size() * sizeof(T));
        buffer_pos += vec.size() * sizeof(T);
    }
    template <class T>
    void serialize(const std::vector<std::vector<T>>& vec, uint8_t *& buffer_pos){
        uint8_t const * const start = buffer_pos;
        for(int i=0; i<vec.size(); i++){
            *reinterpret_cast<uint32_t*>(buffer_pos) = vec[i].size();
            buffer_pos += sizeof(uint32_t);
            memcpy(buffer_pos, vec[i].data(), vec[i].size() * sizeof(T));
            buffer_pos += vec[i].size() * sizeof(T);
        }
    }
    template <class T>
    inline void deserialize(uint8_t const *& buffer_pos, uint32_t size, std::vector<T>& vec){
        vec.clear();
        vec = std::vector<T>(reinterpret_cast<const T*>(buffer_pos), reinterpret_cast<const T*>(buffer_pos) + size);
        buffer_pos += size * sizeof(T);
    }
    template <class T>
    void deserialize(uint8_t const *& buffer_pos, uint32_t num_levels, std::vector<std::vector<T>>& vec){
        vec.clear();
        for(int i=0; i<num_levels; i++){
            uint32_t num = *reinterpret_cast<const uint32_t*>(buffer_pos);
            buffer_pos += sizeof(uint32_t);
            std::vector<T> level_vec = std::vector<T>(reinterpret_cast<const T *>(buffer_pos), reinterpret_cast<const T *>(buffer_pos) + num);
            vec.push_back(level_vec);
            buffer_pos += num * sizeof(T);
        }
    }

    // print vector
    template <class T>
    void print_vec(const std::vector<T>& vec){
        for(int i=0; i<vec.size(); i++){
            std::cout << vec[i] << " ";
        }
        std::cout << std::endl;
    }
    // print nested vector
    template <class T>
    void print_vec(const std::string& name, const std::vector<std::vector<T>>& vec){
        std::cout << name << std::endl;
        for(int i=0; i<vec.size(); i++){
            print_vec(vec[i]);
        }
        std::cout << std::endl;
    }

    // Pack the lowest bit_count bits of each value into a compact byte stream.
    // Stream layout (identical to Jiajun_save_fixed_length_bits from ompSZp_typemanager):
    //   first, for each value, the upper bit_count/8*8 bits (i.e. value >> (bit_count%8)) as little-endian bytes;
    //   then the remaining (bit_count%8)-bit fields of all values, packed as an MSB-first bit stream.
    // Returns the number of bytes written; result must hold at least
    // (bit_count/8)*n + ceil((bit_count%8)*n/8) bytes.
    inline size_t save_fixed_length_bits(unsigned int const * values, size_t n, unsigned char * result, unsigned int bit_count){
        const unsigned int byte_count = bit_count / 8;
        const unsigned int remainder_bit = bit_count % 8;
        const size_t byte_offset = byte_count * n;
        size_t byte_length = byte_offset;
        if(remainder_bit > 0 && n > 0) byte_length += (remainder_bit * n - 1) / 8 + 1;
        // full bytes of each value, little-endian
        if(byte_count > 0){
            unsigned char * pos = result;
            for(size_t i=0; i<n; i++){
                unsigned int tmp = values[i] >> remainder_bit;
                for(unsigned int j=0; j<byte_count; j++){
                    *(pos ++) = (unsigned char) tmp;
                    tmp >>= 8;
                }
            }
        }
        // low remainder_bit bits of each value, MSB-first bit stream
        if(remainder_bit > 0){
            unsigned char * pos = result + byte_offset;
            unsigned char current = 0;
            int bits_filled = 0;
            for(size_t i=0; i<n; i++){
                unsigned int value = values[i] & ((1u << remainder_bit) - 1);
                int bits_left = remainder_bit;
                while(bits_left > 0){
                    int space = 8 - bits_filled;
                    if(bits_left <= space){
                        current |= (unsigned char)(value << (space - bits_left));
                        bits_filled += bits_left;
                        bits_left = 0;
                    }
                    else{
                        current |= (unsigned char)(value >> (bits_left - space));
                        bits_left -= space;
                        value &= (1u << bits_left) - 1;
                        bits_filled = 8;
                    }
                    if(bits_filled == 8){
                        *(pos ++) = current;
                        current = 0;
                        bits_filled = 0;
                    }
                }
            }
            if(bits_filled > 0) *(pos ++) = current;
        }
        return byte_length;
    }

    // Inverse of save_fixed_length_bits: unpack n bit_count-bit values from the byte stream.
    // Returns the number of bytes consumed.
    inline size_t extract_fixed_length_bits(unsigned char const * result, size_t n, unsigned int * values, unsigned int bit_count){
        const unsigned int byte_count = bit_count / 8;
        const unsigned int remainder_bit = bit_count % 8;
        const size_t byte_offset = byte_count * n;
        size_t byte_length = byte_offset;
        if(remainder_bit > 0 && n > 0) byte_length += (remainder_bit * n - 1) / 8 + 1;
        memset(values, 0, n * sizeof(unsigned int));
        // low remainder_bit bits of each value, MSB-first bit stream
        if(remainder_bit > 0){
            unsigned char const * pos = result + byte_offset;
            unsigned char current = 0;
            int bits_available = 0;
            for(size_t i=0; i<n; i++){
                unsigned int value = 0;
                int bits_needed = remainder_bit;
                while(bits_needed > 0){
                    if(bits_available == 0){
                        current = *(pos ++);
                        bits_available = 8;
                    }
                    int take = (bits_available < bits_needed) ? bits_available : bits_needed;
                    value = (value << take) | ((current >> (bits_available - take)) & ((1u << take) - 1));
                    bits_available -= take;
                    bits_needed -= take;
                }
                values[i] = value;
            }
        }
        // full bytes of each value, little-endian
        if(byte_count > 0){
            unsigned char const * pos = result;
            for(size_t i=0; i<n; i++){
                unsigned int tmp = 0;
                for(unsigned int j=0; j<byte_count; j++){
                    tmp |= ((unsigned int) *(pos ++)) << (8 * j);
                }
                values[i] |= tmp << remainder_bit;
            }
        }
        return byte_length;
    }

    class Timer{
    public:
        void start(){
            err = clock_gettime(CLOCK_REALTIME, &start_time);
        }
        void end(){
            err = clock_gettime(CLOCK_REALTIME, &end_time);
            total_time += (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/(double)1000000000;
        }
        double get(){
            double time = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/(double)1000000000;
            clear();
            return time;
        }
        void clear(){
            total_time = 0;
        }
        void print(std::string s){
            std::cout << s << " time: " << total_time << "s" << std::endl;
            clear();
        }
    private:
        int err = 0;
        double total_time = 0;
        struct timespec start_time, end_time;
    };

}
#endif
