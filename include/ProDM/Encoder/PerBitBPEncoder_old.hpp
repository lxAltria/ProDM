#ifndef _MDR_PERBIT_BP_ENCODER_OLD_HPP
#define _MDR_PERBIT_BP_ENCODER_OLD_HPP

#include <cstdlib>

#include <cstring>

#include <vector>

#include <iostream>

#include "BitplaneEncoderInterface.hpp"
#include <bitset>
namespace MDR {
    class BitEncoder_old{
    public:
        BitEncoder_old(uint64_t * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        void encode(uint64_t b){
            buffer += b << position;
            position ++;
            if(position == 64){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        void flush(){
            if(position){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        uint32_t size(){
            return (stream_pos - stream_begin);
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t * stream_pos = NULL;
        uint64_t * stream_begin = NULL;
    };

    class BitDecoder_old{
    public:
        BitDecoder_old(uint64_t const * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        uint8_t decode(){
            if(position == 0){
                buffer = *(stream_pos ++);
                position = 64;
            }
            uint8_t b = buffer & 1u;
            buffer >>= 1;
            position --;
            return b;
        }
        uint32_t size(){
            return (stream_pos - stream_begin);
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t const * stream_pos = NULL;
        uint64_t const * stream_begin = NULL;
    };

    #define PER_BIT_BLOCK_SIZE 1
    // per bit bitplane encoder that encodes data by bit using T_stream type buffer
    template<class T_data, class T_stream>
    class PerBitBPEncoder_old : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        PerBitBPEncoder_old(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder_old: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder_old: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder_old: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder_old: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            // determine block size based on bitplane integer type
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
            }
            std::vector<BitEncoder_old> encoders;
            for(int i=0; i<streams.size(); i++){
                encoders.push_back(BitEncoder_old(reinterpret_cast<uint64_t*>(streams[i])));
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                T_stream sign_bitplane = 0;
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    // compute level errors
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0 && n > 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    // compute level errors
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            for(int i=0; i<num_bitplanes; i++){
                encoders[i].flush();
                stream_sizes[i] = encoders[i].size() * sizeof(uint64_t);
            }
            return streams;
        }

        // std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) {
        //     assert(num_bitplanes > 0);
        //     // determine block size based on bitplane integer type
        //     const int32_t block_size = PER_BIT_BLOCK_SIZE;
        //     stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
        //     std::vector<uint64_t> buffer(num_bitplanes, 0);
        //     std::vector<uint8_t> position(num_bitplanes, 0);
        //     // define fixed point type
        //     using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
        //     std::vector<uint8_t *> streams;
        //     for(int i=0; i<num_bitplanes; i++){
        //         streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
        //     }
        //     std::vector<uint64_t *> stream_begin;
        //     std::vector<uint64_t *> stream_pos;
        //     for(int i=0; i<streams.size(); i++){
        //         stream_begin.push_back(reinterpret_cast<uint64_t*>(streams[i]));
        //         stream_pos.push_back(stream_begin.back());
        //     }
        //     T_data const * data_pos = data;
        //     for(int i=0; i<n - block_size; i+=block_size){
        //         T_stream sign_bitplane = 0;
        //         for(int j=0; j<block_size; j++){
        //             T_data cur_data = *(data_pos++);
        //             T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
        //             bool sign = cur_data < 0;
        //             int64_t fix_point = (int64_t) shifted_data;
        //             T_fp fp_data = sign ? -fix_point : +fix_point;
        //             // compute level errors
        //             bool first_bit = true;
        //             for(int k=num_bitplanes - 1; k>=0; k--){
        //                 uint8_t index = num_bitplanes - 1 - k;
        //                 uint8_t bit = (fp_data >> k) & 1u;
        //                 encode(bit, buffer[index], position[index], stream_pos[index]);
        //                 if(bit && first_bit){
        //                     encode(sign, buffer[index], position[index], stream_pos[index]);
        //                     first_bit = false;
        //                 }
        //             }                    
        //         }
        //     }
        //     // leftover
        //     {
        //         int rest_size = n % block_size;
        //         if(rest_size == 0 && n > 0) rest_size = block_size;
        //         for(int j=0; j<rest_size; j++){
        //             T_data cur_data = *(data_pos++);
        //             T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
        //             bool sign = cur_data < 0;
        //             int64_t fix_point = (int64_t) shifted_data;
        //             T_fp fp_data = sign ? -fix_point : +fix_point;
        //             // compute level errors
        //             bool first_bit = true;
        //             for(int k=num_bitplanes - 1; k>=0; k--){
        //                 uint8_t index = num_bitplanes - 1 - k;
        //                 uint8_t bit = (fp_data >> k) & 1u;
        //                 encode(bit, buffer[index], position[index], stream_pos[index]);
        //                 if(bit && first_bit){
        //                     encode(sign, buffer[index], position[index], stream_pos[index]);
        //                     first_bit = false;
        //                 }
        //             }                    
        //         }
        //     }
        //     for(int i=0; i<num_bitplanes; i++){
        //         flush(buffer[i], position[i], stream_pos[i]);
        //         stream_sizes[i] = size(stream_pos[i], stream_begin[i]) * sizeof(uint64_t);
        //     }
        //     return streams;
        // }

        // only differs in error collection
        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            assert(num_bitplanes > 0);
            // determine block size based on bitplane integer type
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
            }
            std::vector<BitEncoder_old> encoders;
            for(int i=0; i<streams.size(); i++){
                encoders.push_back(BitEncoder_old(reinterpret_cast<uint64_t*>(streams[i])));
            }
            // init level errors
            level_errors.clear();
            level_errors.resize(num_bitplanes + 1);
            for(int i=0; i<level_errors.size(); i++){
                level_errors[i] = 0;
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                T_stream sign_bitplane = 0;
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    // compute level errors
                    collect_level_errors(level_errors, fabs(shifted_data), num_bitplanes);
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0 && n > 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    // compute level errors
                    collect_level_errors(level_errors, fabs(shifted_data), num_bitplanes);
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            for(int i=0; i<num_bitplanes; i++){
                encoders[i].flush();
                stream_sizes[i] = encoders[i].size() * sizeof(uint64_t);
            }
            // translate level errors
            for(int i=0; i<level_errors.size(); i++){
                level_errors[i] = ldexp(level_errors[i], 2*(- num_bitplanes + exp));
            }
            return streams;
        }

        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) {
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                return data;
            }
            std::vector<BitDecoder_old> decoders;
            for(int i=0; i<streams.size(); i++){
                decoders.push_back(BitDecoder_old(reinterpret_cast<uint64_t const*>(streams[i])));
                decoders[i].size();
            }
            // decode
            T_data * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                for(int j=0; j<block_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool first_bit = true;
                    bool sign = false;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = decoders[index].decode();
                        fp_data += bit << k;
                        if(bit && first_bit){
                            // decode sign
                            sign = decoders[index].decode();
                            first_bit = false;
                        }
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - num_bitplanes + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0 && n > 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool first_bit = true;
                    bool sign = false;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = decoders[index].decode();
                        fp_data += bit << k;
                        if(bit && first_bit){
                            // decode sign
                            sign = decoders[index].decode();
                            first_bit = false;
                        }
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - num_bitplanes + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            return data;
        }

        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            return progressive_decode_into(data, streams, n, exp, starting_bitplane, num_bitplanes, level);
        }
        // decode directly into a caller-provided buffer of size n (avoids an internal malloc + external copy)
        T_data * progressive_decode_into(T_data * data, const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            if(num_bitplanes == 0){
                if(level_signs.size() == level){
                    level_signs.push_back(std::vector<bool>(n, false));
                    sign_flags.push_back(std::vector<bool>(n, false));
                }
                memset(data, 0, n * sizeof(T_data));
                return data;
            }
            std::vector<BitDecoder_old> decoders;
            for(int i=0; i<streams.size(); i++){
                decoders.push_back(BitDecoder_old(reinterpret_cast<uint64_t const*>(streams[i])));
                decoders[i].size();
            }
            if(level_signs.size() == level){
                level_signs.push_back(std::vector<bool>(n, false));
                sign_flags.push_back(std::vector<bool>(n, false));
            }
            // std::cout << "level_signs.size() = " << level_signs.size() << ", level = " << level << std::endl;
            std::vector<bool>& signs = level_signs[level];
            std::vector<bool>& flags = sign_flags[level];
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            // decode
            T_data * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                for(int j=0; j<block_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool sign = false;
                    if(flags[i + j]){
                        // sign recorded
                        sign = signs[i + j];
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            T_fp bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        // decode sign if possible
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            T_fp bit = decoders[index].decode();
                            fp_data += bit << k;
                            if(bit && first_bit){
                                // decode sign
                                sign = decoders[index].decode();
                                first_bit = false;
                                flags[i + j] = true;
                            }
                        }
                        signs[i + j] = sign;
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0 && n > 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool sign = false;
                    if(flags[n - rest_size + j]){
                        sign = signs[n - rest_size + j];
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            T_fp bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            T_fp bit = decoders[index].decode();
                            fp_data += bit << k;
                            if(bit && first_bit){
                                // decode sign
                                sign = decoders[index].decode();
                                first_bit = false;
                                flags[n - rest_size + j] = true;
                            }
                        }
                        signs[n - rest_size + j] = sign;
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            return data;
        }

        // T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
        //     const int32_t block_size = PER_BIT_BLOCK_SIZE;
        //     // define fixed point type
        //     using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
        //     T_data * data = (T_data *) malloc(n * sizeof(T_data));
        //     if(num_bitplanes == 0){
        //         if(level_signs.size() == level){
        //             level_signs.push_back(std::vector<bool>(n, false));
        //             sign_flags.push_back(std::vector<bool>(n, false));
        //         }
        //         memset(data, 0, n * sizeof(T_data));
        //         return data;
        //     }
        //     std::vector<uint64_t> buffer(num_bitplanes, 0);
        //     std::vector<uint8_t> position(num_bitplanes, 0);
        //     // std::vector<BitDecoder_old> decoders;
        //     std::vector<uint64_t const *> stream_pos;
        //     std::vector<uint64_t const *> stream_begin;
        //     for(int i=0; i<streams.size(); i++){
        //         stream_begin.push_back(reinterpret_cast<uint64_t const*>(streams[i]));
        //         stream_pos.push_back(stream_begin.back());
        //         // decoders.push_back(BitDecoder_old(reinterpret_cast<uint64_t const*>(streams[i])));
        //         // decoders[i].size();
        //         size(stream_pos[i], stream_begin[i]);
        //     }
        //     if(level_signs.size() == level){
        //         level_signs.push_back(std::vector<bool>(n, false));
        //         sign_flags.push_back(std::vector<bool>(n, false));
        //     }
        //     // std::cout << "level_signs.size() = " << level_signs.size() << ", level = " << level << std::endl;
        //     std::vector<bool>& signs = level_signs[level];
        //     std::vector<bool>& flags = sign_flags[level];
        //     const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
        //     // decode
        //     T_data * data_pos = data;
        //     for(int i=0; i<n - block_size; i+=block_size){
        //         for(int j=0; j<block_size; j++){
        //             T_fp fp_data = 0;
        //             // decode each bit of the data for each level component
        //             bool sign = false;
        //             if(flags[i + j]){
        //                 // sign recorded
        //                 sign = signs[i + j];
        //                 for(int k=num_bitplanes - 1; k>=0; k--){
        //                     uint8_t index = num_bitplanes - 1 - k;
        //                     T_fp bit = decode(buffer[index], position[index], stream_pos[index]);
        //                     fp_data += bit << k;
        //                 }
        //             }
        //             else{
        //                 // decode sign if possible
        //                 bool first_bit = true;
        //                 for(int k=num_bitplanes - 1; k>=0; k--){
        //                     uint8_t index = num_bitplanes - 1 - k;
        //                     T_fp bit = decode(buffer[index], position[index], stream_pos[index]);
        //                     fp_data += bit << k;
        //                     if(bit && first_bit){
        //                         // decode sign
        //                         sign = decode(buffer[index], position[index], stream_pos[index]);
        //                         first_bit = false;
        //                         flags[i + j] = true;
        //                     }
        //                 }
        //                 signs[i + j] = sign;
        //             }
        //             T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
        //             *(data_pos++) = sign ? -cur_data : cur_data;
        //         }
        //     }
        //     // leftover
        //     {
        //         int rest_size = n % block_size;
        //         if(rest_size == 0 && n > 0) rest_size = block_size;
        //         for(int j=0; j<rest_size; j++){
        //             T_fp fp_data = 0;
        //             // decode each bit of the data for each level component
        //             bool sign = false;
        //             if(flags[n - rest_size + j]){
        //                 sign = signs[n - rest_size + j];
        //                 for(int k=num_bitplanes - 1; k>=0; k--){
        //                     uint8_t index = num_bitplanes - 1 - k;
        //                     T_fp bit = decode(buffer[index], position[index], stream_pos[index]);
        //                     fp_data += bit << k;
        //                 }
        //             }
        //             else{
        //                 bool first_bit = true;
        //                 for(int k=num_bitplanes - 1; k>=0; k--){
        //                     uint8_t index = num_bitplanes - 1 - k;
        //                     T_fp bit = decode(buffer[index], position[index], stream_pos[index]);
        //                     fp_data += bit << k;
        //                     if(bit && first_bit){
        //                         // decode sign
        //                         sign = decode(buffer[index], position[index], stream_pos[index]);
        //                         first_bit = false;
        //                         flags[n - rest_size + j] = true;
        //                     }
        //                 }
        //                 signs[n - rest_size + j] = sign;
        //             }
        //             T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
        //             *(data_pos++) = sign ? -cur_data : cur_data;
        //         }
        //     }
        //     return data;
        // }

        void print() const {
            std::cout << "Per-bit bitplane encoder" << std::endl;
        }
    private:
        inline void encode(uint64_t b, uint64_t & buffer, uint8_t & position, uint64_t * & stream_pos){
            buffer += b << position;
            position ++;
            if(position == 64){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        inline void flush(uint64_t & buffer, uint8_t & position, uint64_t * & stream_pos){
            if(position){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        inline uint32_t size(uint64_t * & stream_pos, uint64_t * & stream_begin){
            return (stream_pos - stream_begin);
        }
        inline uint32_t size(uint64_t const * & stream_pos, uint64_t const * & stream_begin){
            return (stream_pos - stream_begin);
        }
        inline uint8_t decode(uint64_t & buffer, uint8_t & position, uint64_t const * & stream_pos){
            if(position == 0){
                buffer = *(stream_pos ++);
                position = 64;
            }
            uint8_t b = buffer & 1u;
            buffer >>= 1;
            position --;
            return b;
        }
        inline void collect_level_errors(std::vector<double>& level_errors, float data, int num_bitplanes) const {
            uint32_t fp_data = (uint32_t) data;
            double mantissa = data - (uint32_t) data;
            level_errors[num_bitplanes] += mantissa * mantissa;
            for(int k=1; k<num_bitplanes; k++){
                uint32_t mask = (1 << k) - 1;
                double diff = (double) (fp_data & mask) + mantissa;
                level_errors[num_bitplanes - k] += diff * diff;
            }
            level_errors[0] += data * data;
        }
        std::vector<std::vector<bool>> level_signs;
        std::vector<std::vector<bool>> sign_flags;
    };
}
#endif
