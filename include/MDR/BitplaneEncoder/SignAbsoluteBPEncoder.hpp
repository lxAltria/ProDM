#ifndef _SIGN_ABSOLUTE_BP_ENCODER_HPP
#define _SIGN_ABSOLUTE_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"

namespace MDR {
    template<class T_data, class T_stream>
    class SignAbsoluteBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        SignAbsoluteBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "SignAbsoluteBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "SignAbsoluteBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "SignAbsoluteBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "SignAbsoluteBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            stream_sizes = std::vector<uint32_t>(num_bitplanes + 1, 0);

            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;

            streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            for(int i = 0; i < num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }

            std::vector<T_stream *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }

            std::vector<T_fp> int_data_buffer(block_size);
            T_data const * data_pos = data;
            T_stream prev_sign_word = 0;

            for(int i = 0; i < n - (int)block_size; i += block_size){
                for(uint32_t j = 0; j < block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = (T_fp) fabs(shifted_data);
                }
                encode_block_with_xor(int_data_buffer.data(), data_pos - block_size, block_size, num_bitplanes, streams_pos, prev_sign_word);
            }

            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j = 0; j < rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = (T_fp) fabs(shifted_data);
                }
                encode_block_with_xor(int_data_buffer.data(), data_pos - rest_size, rest_size, num_bitplanes, streams_pos, prev_sign_word);
            }

            for(int i = 0; i < (int)streams.size(); i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            return streams;
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            assert(num_bitplanes > 0);
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            stream_sizes = std::vector<uint32_t>(num_bitplanes + 1, 0);

            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;

            streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            for(int i = 0; i < num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }

            std::vector<T_stream *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }

            level_errors.clear();
            level_errors.resize(num_bitplanes + 1, 0.0);

            std::vector<T_fp> int_data_buffer(block_size);
            T_data const * data_pos = data;
            T_stream prev_sign_word = 0;

            for(int i = 0; i < n - (int)block_size; i += block_size){
                for(uint32_t j = 0; j < block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = (T_fp) fabs(shifted_data);
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, num_bitplanes);
                }
                encode_block_with_xor(int_data_buffer.data(), data_pos - block_size, block_size, num_bitplanes, streams_pos, prev_sign_word);
            }

            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j = 0; j < rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = (T_fp) fabs(shifted_data);
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, num_bitplanes);
                }
                encode_block_with_xor(int_data_buffer.data(), data_pos - rest_size, rest_size, num_bitplanes, streams_pos, prev_sign_word);
            }

            for(int i = 0; i < (int)streams.size(); i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            for(size_t i = 0; i < level_errors.size(); i++){
                level_errors[i] = ldexp(level_errors[i], 2 * (- (int)num_bitplanes + exp));
            }
            return streams;
        }

        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) {
            return progressive_decode(streams, n, exp, 0, num_bitplanes, 0);
        }

        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;

            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                return data;
            }

            std::vector<T_stream const *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream const *>(streams[i]);
            }

            std::vector<T_fp> int_data_buffer(block_size);
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            const T_data scale = ldexp((T_data)1.0, -(int)ending_bitplane + exp);

            T_data * data_pos = data;
            T_stream prev_sign_word = 0;

            for(int i = 0; i < n - (int)block_size; i += (int)block_size){
                T_stream sign_word_encoded = *(streams_pos[0]++);
                T_stream sign_word = sign_word_encoded ^ prev_sign_word;
                prev_sign_word = sign_word;

                memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));

                for(int k = num_bitplanes - 1; k >= 0; k--){
                    T_stream encoded = *(streams_pos[k + 1]++);

                    // Recover from difference encoding
                    T_stream recovered = encoded & 1u;  // first bit unchanged
                    for(uint32_t j = 1; j < block_size; j++){
                        T_stream curr_encoded = (encoded >> j) & 1u;
                        T_stream prev_recovered = (recovered >> (j - 1)) & 1u;
                        recovered |= (curr_encoded ^ prev_recovered) << j;
                    }

                    for(uint32_t j = 0; j < block_size; j++){
                        int_data_buffer[j] |= (T_fp)((recovered >> j) & 1u) << k;
                    }
                }

                for(uint32_t j = 0; j < block_size; j++){
                    bool sign = (sign_word >> j) & 1u;
                    T_data magnitude = scale * (T_data)int_data_buffer[j];
                    *(data_pos++) = sign ? -magnitude : magnitude;
                }
            }

            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                T_stream sign_word_encoded = *(streams_pos[0]++);
                T_stream sign_word = sign_word_encoded ^ prev_sign_word;

                memset(int_data_buffer.data(), 0, rest_size * sizeof(T_fp));

                for(int k = num_bitplanes - 1; k >= 0; k--){
                    T_stream encoded = *(streams_pos[k + 1]++);

                    // Recover from difference encoding
                    T_stream recovered = encoded & 1u;  // first bit unchanged
                    for(int j = 1; j < rest_size; j++){
                        T_stream curr_encoded = (encoded >> j) & 1u;
                        T_stream prev_recovered = (recovered >> (j - 1)) & 1u;
                        recovered |= (curr_encoded ^ prev_recovered) << j;
                    }

                    for(int j = 0; j < rest_size; j++){
                        int_data_buffer[j] |= (T_fp)((recovered >> j) & 1u) << k;
                    }
                }

                for(int j = 0; j < rest_size; j++){
                    bool sign = (sign_word >> j) & 1u;
                    T_data magnitude = scale * (T_data)int_data_buffer[j];
                    *(data_pos++) = sign ? -magnitude : magnitude;
                }
            }

            return data;
        }

        void print() const {
            std::cout << "Sign-Absolute bitplane encoder (XOR + separated signs)" << std::endl;
        }

    private:
        template<class T>
        uint32_t block_size_based_on_bitplane_int_type() const {
            if(std::is_same<T, uint64_t>::value) return 64;
            if(std::is_same<T, uint32_t>::value) return 32;
            if(std::is_same<T, uint16_t>::value) return 16;
            if(std::is_same<T, uint8_t>::value) return 8;
            std::cerr << "Integer type not supported." << std::endl;
            exit(0);
        }

        template <class T_int>
        inline void encode_block_with_xor(T_int const * magnitude, T_data const * original, size_t n, uint8_t num_bitplanes, std::vector<T_stream *>& streams_pos, T_stream& prev_sign_word) const {
            T_stream sign_word = 0;
            for(size_t i = 0; i < n; i++){
                bool sign = original[i] < 0;
                sign_word |= (T_stream)sign << i;
            }
            T_stream sign_word_encoded = sign_word ^ prev_sign_word;
            *(streams_pos[0]++) = sign_word_encoded;
            prev_sign_word = sign_word;

            // Intra-bp XOR: within same bitplane, adjacent elements difference
            for(int k = num_bitplanes - 1; k >= 0; k--){
                T_stream bitplane_value = 0;
                for(size_t i = 0; i < n; i++){
                    T_stream bit_k = (T_stream)((magnitude[i] >> k) & 1u);
                    bitplane_value |= bit_k << i;
                }

                // Difference encoding: bit[i] XOR bit[i-1]
                T_stream encoded = bitplane_value & 1u;  // first bit unchanged
                for(size_t i = 1; i < n; i++){
                    T_stream bit_curr = (bitplane_value >> i) & 1u;
                    T_stream bit_prev = (bitplane_value >> (i - 1)) & 1u;
                    encoded |= (bit_curr ^ bit_prev) << i;
                }

                *(streams_pos[k + 1]++) = encoded;
            }
        }

        inline void collect_level_errors(std::vector<double>& level_errors, uint32_t magnitude, float data, int num_bitplanes) const {
            float abs_data = fabs(data);
            double mantissa = abs_data - (uint32_t) abs_data;
            level_errors[num_bitplanes] += mantissa * mantissa;
            for(int k = 1; k < num_bitplanes; k++){
                uint32_t mask = (1u << k) - 1;
                double diff = (double)(magnitude & mask) + mantissa;
                level_errors[num_bitplanes - k] += diff * diff;
            }
            level_errors[0] += abs_data * abs_data;
        }
    };
}
#endif
