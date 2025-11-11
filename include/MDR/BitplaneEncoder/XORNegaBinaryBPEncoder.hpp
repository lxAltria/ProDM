#ifndef _MDR_XORNEGABINARY_BP_ENCODER_HPP
#define _MDR_XORNEGABINARY_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include <unordered_set>

namespace MDR {
    // general bitplane encoder that encodes data by block using T_stream type buffer
    template<class T_data, class T_stream>
    class XORNegaBinaryBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        XORNegaBinaryBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "XORNegaBinaryBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "XORNegaBinaryBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "XORNegaBinaryBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "XORNegaBinaryBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            // leave room for negabinary format
            exp += 2;
            // determine block size based on bitplane integer type
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            std::vector<uint8_t> starting_bitplanes = std::vector<uint8_t>((n - 1)/block_size + 1, 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            std::vector<T_stream *> streams_pos(streams.size());
            for(int i=0; i<streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = binary2negabinary((T_fps) shifted_data);
                }
                encode_block(int_data_buffer.data(), block_size, num_bitplanes, streams_pos);
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = binary2negabinary((T_fps) shifted_data);
                }
                encode_block(int_data_buffer.data(), rest_size, num_bitplanes, streams_pos);
            }
            for(int i=0; i<num_bitplanes; i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            return streams;
        }

        // only differs in error collection
        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            assert(num_bitplanes > 0);
            // leave room for negabinary format
            exp += 2;
            // determine block size based on bitplane integer type
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            std::vector<uint8_t> starting_bitplanes = std::vector<uint8_t>((n - 1)/block_size + 1, 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            std::vector<T_stream *> streams_pos(streams.size());
            for(int i=0; i<streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }
            // init level errors
            level_errors.clear();
            level_errors.resize(num_bitplanes + 1);
            for(int i=0; i<level_errors.size(); i++){
                level_errors[i] = 0;
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_fps signed_int_data = (T_fps) shifted_data;
                    int_data_buffer[j] = binary2negabinary(signed_int_data);
                    // compute level errors
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, shifted_data - signed_int_data, num_bitplanes);
                }
                encode_block(int_data_buffer.data(), block_size, num_bitplanes, streams_pos);
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_fps signed_int_data = (T_fps) shifted_data;
                    int_data_buffer[j] = binary2negabinary(signed_int_data);
                    // compute level errors
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, shifted_data - signed_int_data, num_bitplanes);
                }
                encode_block(int_data_buffer.data(), rest_size, num_bitplanes, streams_pos);
            }
            for(int i=0; i<num_bitplanes; i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            // translate level errors
            for(int i=0; i<level_errors.size(); i++){
                level_errors[i] = ldexp(level_errors[i], 2*(- num_bitplanes + exp));
            }
            return streams;
        }

        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) {
            return progressive_decode(streams, n, exp, 0, num_bitplanes, streams.size());
        }

        // decode the data and record necessary information for progressiveness
        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            bool first_time_decode = (level_table.count(level)) ? false : true;
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                return data;
            }
            // leave room for negabinary format
            exp += 2;
            // define fixed point type
            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<T_stream const *> streams_pos(streams.size());
            for(int i=0; i<streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream const *>(streams[i]);
            }
            if(first_time_decode){
                if(last_2_data.size() <= level){
                    last_2_data.resize(level + 1);
                    last_1_data.resize(level + 1);
                }
                for(int i=0; i<n - block_size; i+=block_size){
                    last_2_data[level].push_back(std::vector<unsigned char>(block_size));
                    last_1_data[level].push_back(std::vector<unsigned char>(block_size));
                }
                {
                    int rest_size = n % block_size;
                    if(!rest_size) rest_size = block_size;
                    last_2_data[level].push_back(std::vector<unsigned char>(rest_size));
                    last_1_data[level].push_back(std::vector<unsigned char>(rest_size));
                }
                level_table.insert(level);
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            // decode
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            T_data * data_pos = data;
            // std::cout << "ending_bitplane = " << +ending_bitplane << std::endl;
            if(ending_bitplane % 2 == 0){
                int block_idx = 0;
                for(int i=0; i<n - block_size; i+=block_size){
                    memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
                    decode_block(streams_pos, block_size, num_bitplanes, int_data_buffer.data(), level, block_idx, first_time_decode);
                    for(int j=0; j<block_size; j++){
                        *(data_pos++) = ldexp((T_data) negabinary2binary(int_data_buffer[j]), - ending_bitplane + exp);
                    }
                    block_idx++;
                }
                // leftover
                {
                    int rest_size = n % block_size;
                    if(rest_size == 0) rest_size = block_size;
                    memset(int_data_buffer.data(), 0, rest_size * sizeof(T_fp));
                    decode_block(streams_pos, rest_size, num_bitplanes, int_data_buffer.data(), level, block_idx, first_time_decode);
                    for(int j=0; j<rest_size; j++){
                        *(data_pos++) = ldexp((T_data) negabinary2binary(int_data_buffer[j]), - ending_bitplane + exp);
                    }
                }                
            }
            else{
                int block_idx = 0;
                for(int i=0; i<n - block_size; i+=block_size){
                    memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
                    decode_block(streams_pos, block_size, num_bitplanes, int_data_buffer.data(), level, block_idx, first_time_decode);
                    for(int j=0; j<block_size; j++){
                        *(data_pos++) = - ldexp((T_data) negabinary2binary(int_data_buffer[j]), - ending_bitplane + exp);
                    }
                    block_idx++;
                }
                // leftover
                {
                    int rest_size = n % block_size;
                    if(rest_size == 0) rest_size = block_size;
                    memset(int_data_buffer.data(), 0, rest_size * sizeof(T_fp));
                    decode_block(streams_pos, rest_size, num_bitplanes, int_data_buffer.data(), level, block_idx, first_time_decode);
                    for(int j=0; j<rest_size; j++){
                        *(data_pos++) = - ldexp((T_data) negabinary2binary(int_data_buffer[j]), - ending_bitplane + exp);
                    }
                }                
            }
            return data;
        }

        void print() const {
            std::cout << "XORNegaBinary bitplane encoder" << std::endl;
        }
    private:
        template<class T>
        uint32_t block_size_based_on_bitplane_int_type() const {
            uint32_t block_size = 0;
            if(std::is_same<T, uint64_t>::value){
                block_size = 64;
            }
            else if(std::is_same<T, uint32_t>::value){
                block_size = 32;
            }
            else if(std::is_same<T, uint16_t>::value){
                block_size = 16;
            }
            else if(std::is_same<T, uint8_t>::value){
                block_size = 8;
            }
            else{
                std::cerr << "Integer type not supported." << std::endl;
                exit(0);
            }
            return block_size;
        }
        inline uint64_t binary2negabinary(const int64_t x) const {
            return (x + (uint64_t)0xaaaaaaaaaaaaaaaaull) ^ (uint64_t)0xaaaaaaaaaaaaaaaaull;
        }
        inline uint32_t binary2negabinary(const int32_t x) const {
            return (x + (uint32_t)0xaaaaaaaau) ^ (uint32_t)0xaaaaaaaau;
        }
        inline int64_t negabinary2binary(const uint64_t x) const {
            return (x ^0xaaaaaaaaaaaaaaaaull) - 0xaaaaaaaaaaaaaaaaull;
        }
        inline int32_t negabinary2binary(const uint32_t x) const {
            return (x ^0xaaaaaaaau) - 0xaaaaaaaau;
        }
        inline void collect_level_errors(std::vector<double>& level_errors, uint32_t negabinary_data, float data, float mantissa, int num_bitplanes) const {
            level_errors[num_bitplanes] += mantissa * mantissa;
            for(int k=1; k<num_bitplanes; k++){
                uint32_t mask = (1 << k) - 1;
                double diff = (double) negabinary2binary(negabinary_data & mask) + mantissa;
                level_errors[num_bitplanes - k] += diff * diff;
            }
            level_errors[0] += data * data;
        }
        template <class T_int>
        inline void encode_block(T_int const * data, size_t n, uint8_t num_bitplanes, std::vector<T_stream *>& streams_pos) const {
            for(int k=num_bitplanes - 1; k>=0; k--){
                T_stream bitplane_value = 0;
                T_stream bitplane_index = num_bitplanes - 1 - k;
                if(k < num_bitplanes - 2){
                    for (int i=0; i<n; i++){
                        T_stream last_2_bit = (T_stream)((data[i] >> (k + 2)) & 1u);
                        T_stream last_1_bit = (T_stream)((data[i] >> (k + 1)) & 1u);
                        T_stream current_bit = (T_stream)((data[i] >> k) & 1u);
                        bitplane_value += (T_stream)(last_2_bit ^ last_1_bit ^ current_bit) << i;
                    }
                }
                else{
                    for (int i=0; i<n; i++){
                        bitplane_value += (T_stream)((data[i] >> k) & 1u) << i;
                    }
                }
                *(streams_pos[bitplane_index] ++) = bitplane_value;
            }
        }
        template <class T_int>
        inline void decode_block(std::vector<T_stream const *>& streams_pos, size_t n, uint8_t num_bitplanes, T_int * data, int level, int block_index, bool first_time_decode) const {
            // First time decode 1 bitplane has not been considered. But it is rare.
            for(int k=num_bitplanes - 1; k>=0; k--){
                T_stream bitplane_index = num_bitplanes - 1 - k;
                T_stream bitplane_value = *(streams_pos[bitplane_index] ++);
                if(first_time_decode && k >= num_bitplanes - 2){
                    for (int i=0; i<n; i++){
                        data[i] += ((bitplane_value >> i) & 1u) << k;
                    }
                }
                else {
                    for (int i=0; i<n; i++){
                        T_stream last_2_bit = (T_stream)(last_2_data[level][block_index][i]);
                        T_stream last_1_bit = (T_stream)(last_1_data[level][block_index][i]);
                        T_stream temp = last_2_bit ^ last_1_bit;
                        T_stream current_bit = ((bitplane_value >> i) & 1u);
                        T_stream recovered_bit = (!current_bit) ? temp : (1-temp);
                        data[i] += recovered_bit << k;
                    }
                }
                for (int i=0; i<n; i++){
                    last_2_data[level][block_index][i] = last_1_data[level][block_index][i];
                    last_1_data[level][block_index][i] = static_cast<unsigned char>((data[i] >> k) & 1u);
                }
            }
        }

        mutable std::unordered_set<int> level_table;
        mutable std::vector<std::vector<std::vector<unsigned char>>> last_2_data; // [level][block_idx][i]
        mutable std::vector<std::vector<std::vector<unsigned char>>> last_1_data;
    };
}
#endif
