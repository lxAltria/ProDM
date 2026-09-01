#ifndef _MDR_PERBIT_BP_ENCODER_HPP
#define _MDR_PERBIT_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include <bitset>
namespace MDR {
    class BitEncoder{
    public:
        BitEncoder(uint64_t * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        inline void encode(uint64_t b){
            buffer += b << position;
            position ++;
            if(position == 64){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        // Write exactly 64 data bits at once, bypassing per-bit branching
        inline void write_aligned_64(uint64_t word){
            if(position == 0){
                *(stream_pos++) = word;
            } else {
                buffer |= (word << position);
                *(stream_pos++) = buffer;
                buffer = word >> (64 - position);
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

    class BitDecoder{
    public:
        BitDecoder(uint64_t const * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        inline uint8_t decode(){
            if(position == 0){
                buffer = *(stream_pos ++);
                position = 64;
            }
            uint8_t b = buffer & 1u;
            buffer >>= 1;
            position --;
            return b;
        }
        // Read exactly 64 data bits at once, bypassing per-bit branching
        inline uint64_t read_aligned_64(){
            if(position == 0){
                return *(stream_pos++);
            } else {
                uint64_t result = buffer;
                uint64_t next = *(stream_pos++);
                result |= (next << position);
                buffer = next >> (64 - position);
                return result;
            }
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
    class PerBitBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        PerBitBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
            }
            std::vector<BitEncoder> encoders;
            for(int i=0; i<(int)streams.size(); i++){
                encoders.push_back(BitEncoder(reinterpret_cast<uint64_t*>(streams[i])));
            }
            // Tile buffers — sized to fit in L2 cache
            constexpr int32_t TILE = 4096;
            std::vector<T_fp> fp_buf(TILE);
            std::vector<uint8_t> sign_buf(TILE);
            std::vector<uint8_t> lead_bp(TILE);
            for(int32_t t = 0; t < n; t += TILE){
                int32_t tlen = std::min(TILE, n - t);
                // Pass 1: convert to fixed-point, find leading bitplane via CLZ
                for(int32_t i = 0; i < tlen; i++){
                    T_data cur = data[t + i];
                    T_data shifted = ldexp(cur, num_bitplanes - exp);
                    bool s = cur < 0;
                    int64_t fix = (int64_t) shifted;
                    T_fp fp = s ? (T_fp)(-fix) : (T_fp)(fix);
                    fp_buf[i] = fp;
                    sign_buf[i] = s ? 1 : 0;
                    if(fp != 0){
                        int hi = (sizeof(T_fp) == 8) ? 63 - __builtin_clzll((uint64_t)fp)
                                                     : 31 - __builtin_clz((uint32_t)fp);
                        lead_bp[i] = (hi < num_bitplanes) ? (uint8_t)(num_bitplanes - 1 - hi) : num_bitplanes;
                    } else {
                        lead_bp[i] = num_bitplanes;
                    }
                }
                // Pass 2: bitplane-major — one encoder at a time (cache-friendly)
                for(uint8_t bp = 0; bp < num_bitplanes; bp++){
                    int k = num_bitplanes - 1 - bp;
                    int32_t i = 0;
                    // Fast path: batch 64 bits into one word when no sign bits in chunk
                    while(i + 64 <= tlen){
                        bool has_sign = false;
                        for(int32_t j = i; j < i + 64; j++){
                            if(lead_bp[j] == bp){ has_sign = true; break; }
                        }
                        if(!has_sign){
                            uint64_t word = 0;
                            for(int32_t j = 0; j < 64; j++){
                                word |= ((uint64_t)((fp_buf[i+j] >> k) & 1u)) << j;
                            }
                            encoders[bp].write_aligned_64(word);
                        } else {
                            for(int32_t j = i; j < i + 64; j++){
                                uint8_t bit = (fp_buf[j] >> k) & 1u;
                                encoders[bp].encode(bit);
                                if(lead_bp[j] == bp) encoders[bp].encode(sign_buf[j]);
                            }
                        }
                        i += 64;
                    }
                    for(; i < tlen; i++){
                        uint8_t bit = (fp_buf[i] >> k) & 1u;
                        encoders[bp].encode(bit);
                        if(lead_bp[i] == bp) encoders[bp].encode(sign_buf[i]);
                    }
                }
            }
            for(int i=0; i<num_bitplanes; i++){
                encoders[i].flush();
                stream_sizes[i] = encoders[i].size() * sizeof(uint64_t);
            }
            return streams;
        }

        // encode with level error collection
        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            assert(num_bitplanes > 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
            }
            std::vector<BitEncoder> encoders;
            for(int i=0; i<(int)streams.size(); i++){
                encoders.push_back(BitEncoder(reinterpret_cast<uint64_t*>(streams[i])));
            }
            level_errors.clear();
            level_errors.resize(num_bitplanes + 1, 0.0);
            constexpr int32_t TILE = 4096;
            std::vector<T_fp> fp_buf(TILE);
            std::vector<uint8_t> sign_buf(TILE);
            std::vector<uint8_t> lead_bp(TILE);
            for(int32_t t = 0; t < n; t += TILE){
                int32_t tlen = std::min(TILE, n - t);
                for(int32_t i = 0; i < tlen; i++){
                    T_data cur = data[t + i];
                    T_data shifted = ldexp(cur, num_bitplanes - exp);
                    bool s = cur < 0;
                    int64_t fix = (int64_t) shifted;
                    T_fp fp = s ? (T_fp)(-fix) : (T_fp)(fix);
                    fp_buf[i] = fp;
                    sign_buf[i] = s ? 1 : 0;
                    if(fp != 0){
                        int hi = (sizeof(T_fp) == 8) ? 63 - __builtin_clzll((uint64_t)fp)
                                                     : 31 - __builtin_clz((uint32_t)fp);
                        lead_bp[i] = (hi < num_bitplanes) ? (uint8_t)(num_bitplanes - 1 - hi) : num_bitplanes;
                    } else {
                        lead_bp[i] = num_bitplanes;
                    }
                    collect_level_errors(level_errors, fabs(shifted), num_bitplanes);
                }
                for(uint8_t bp = 0; bp < num_bitplanes; bp++){
                    int k = num_bitplanes - 1 - bp;
                    int32_t i = 0;
                    while(i + 64 <= tlen){
                        bool has_sign = false;
                        for(int32_t j = i; j < i + 64; j++){
                            if(lead_bp[j] == bp){ has_sign = true; break; }
                        }
                        if(!has_sign){
                            uint64_t word = 0;
                            for(int32_t j = 0; j < 64; j++){
                                word |= ((uint64_t)((fp_buf[i+j] >> k) & 1u)) << j;
                            }
                            encoders[bp].write_aligned_64(word);
                        } else {
                            for(int32_t j = i; j < i + 64; j++){
                                uint8_t bit = (fp_buf[j] >> k) & 1u;
                                encoders[bp].encode(bit);
                                if(lead_bp[j] == bp) encoders[bp].encode(sign_buf[j]);
                            }
                        }
                        i += 64;
                    }
                    for(; i < tlen; i++){
                        uint8_t bit = (fp_buf[i] >> k) & 1u;
                        encoders[bp].encode(bit);
                        if(lead_bp[i] == bp) encoders[bp].encode(sign_buf[i]);
                    }
                }
            }
            for(int i=0; i<num_bitplanes; i++){
                encoders[i].flush();
                stream_sizes[i] = encoders[i].size() * sizeof(uint64_t);
            }
            for(int i=0; i<(int)level_errors.size(); i++){
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
            std::vector<BitDecoder> decoders;
            for(int i=0; i<streams.size(); i++){
                decoders.push_back(BitDecoder(reinterpret_cast<uint64_t const*>(streams[i])));
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
                if(rest_size == 0) rest_size = block_size;
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
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                if(level_signs.size() == level){
                    level_signs.push_back(std::vector<bool>(n, false));
                    sign_flags.push_back(std::vector<bool>(n, false));
                }
                memset(data, 0, n * sizeof(T_data));
                return data;
            }
            std::vector<BitDecoder> decoders;
            for(int i=0; i<streams.size(); i++){
                decoders.push_back(BitDecoder(reinterpret_cast<uint64_t const*>(streams[i])));
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
                if(rest_size == 0) rest_size = block_size;
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

        void print() const {
            std::cout << "Per-bit bitplane encoder" << std::endl;
        }
    private:
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
