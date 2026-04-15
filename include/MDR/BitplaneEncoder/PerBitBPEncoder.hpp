#ifndef _MDR_PERBIT_BP_ENCODER_HPP
#define _MDR_PERBIT_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include <bitset>
#include <cstring>
#include <cmath>
#include <algorithm>

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

    template<class T_data, class T_stream>
    class PerBitBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        PerBitBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
        }

        // =====================================================================
        // encode (without level errors)
        // =====================================================================
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
            return streams;
        }

        // =====================================================================
        // encode (with level error collection)
        // =====================================================================
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

        // =====================================================================
        // decode (full, non-progressive)
        // =====================================================================
        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) {
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                return data;
            }
            std::vector<BitDecoder> decoders;
            for(size_t i = 0; i < streams.size(); i++){
                decoders.push_back(BitDecoder(reinterpret_cast<uint64_t const*>(streams[i])));
            }
            // Keep element-major order — matches encoder's stream layout
            // Use raw pointer array to avoid vector<BitDecoder> operator[] overhead
            BitDecoder * dec = decoders.data();
            const T_data scale = ldexp((T_data)1.0, -num_bitplanes + exp);
            for(int32_t i = 0; i < n; i++){
                T_fp fp_data = 0;
                bool first_bit = true;
                bool sign = false;
                for(int k = num_bitplanes - 1; k >= 0; k--){
                    uint8_t index = num_bitplanes - 1 - k;
                    uint8_t bit = dec[index].decode();
                    fp_data |= (T_fp)bit << k;
                    if(bit && first_bit){
                        sign = dec[index].decode();
                        first_bit = false;
                    }
                }
                T_data cur_data = scale * (T_data)fp_data;
                data[i] = sign ? -cur_data : cur_data;
            }
            return data;
        }

        // =====================================================================
        // progressive_decode — conservatively optimized
        //
        //  Kept element-major order (must match encoder stream layout).
        //  Changes:
        //   1. vector<bool> → vector<uint8_t> for signs/flags
        //      (avoids bit-packing proxy objects, ~2× faster random access)
        //   2. Raw pointer to decoder array (skip vector operator[])
        //   3. Precompute ldexp scale factor outside the loop
        //   4. Eliminated block_size=1 boilerplate (single flat loop)
        //   5. Use |= instead of += for bit assembly (no carry propagation)
        //   6. Removed dead decoders[i].size() calls
        // =====================================================================
        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp,
                                     uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            T_data * data = (T_data *) malloc(n * sizeof(T_data));

            // Ensure sign/flag vectors exist for this level (uint8_t, not bool)
            if((int)level_signs.size() <= level){
                level_signs.resize(level + 1);
                sign_flags.resize(level + 1);
            }
            if(level_signs[level].empty()){
                level_signs[level].assign(n, 0);
                sign_flags[level].assign(n, 0);
            }

            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                return data;
            }

            std::vector<BitDecoder> decoders;
            decoders.reserve(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                decoders.push_back(BitDecoder(reinterpret_cast<uint64_t const*>(streams[i])));
            }
            // Raw pointer — avoids vector bounds-check / indirection in hot loop
            BitDecoder * dec = decoders.data();

            uint8_t * signs = level_signs[level].data();
            uint8_t * flags = sign_flags[level].data();
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            const T_data scale = ldexp((T_data)1.0, -(int)ending_bitplane + exp);

            // Element-major decode — matches encoder's interleaved stream order
            for(int32_t i = 0; i < n; i++){
                T_fp fp_data = 0;

                if(flags[i]){
                    // Sign already known from a previous progressive round.
                    // Pure data bits — tight loop, no branching on bit value.
                    for(int k = num_bitplanes - 1; k >= 0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        fp_data |= (T_fp)dec[index].decode() << k;
                    }
                    T_data cur_data = scale * (T_data)fp_data;
                    data[i] = signs[i] ? -cur_data : cur_data;
                }
                else{
                    // Sign not yet found — need to watch for first set bit.
                    bool sign = false;
                    bool first_bit = true;
                    for(int k = num_bitplanes - 1; k >= 0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        T_fp bit = dec[index].decode();
                        fp_data |= bit << k;
                        if(bit && first_bit){
                            sign = dec[index].decode();
                            first_bit = false;
                            flags[i] = 1;
                        }
                    }
                    signs[i] = sign ? 1 : 0;
                    T_data cur_data = scale * (T_data)fp_data;
                    data[i] = sign ? -cur_data : cur_data;
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

        // uint8_t instead of bool — no bit-packing, direct byte access
        std::vector<std::vector<uint8_t>> level_signs;
        std::vector<std::vector<uint8_t>> sign_flags;
    };
}
#endif