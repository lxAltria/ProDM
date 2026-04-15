#ifndef _MDR_XORNEGABINARY_BP_ENCODER_HPP
#define _MDR_XORNEGABINARY_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include <cstring>
#include <cmath>
#include <algorithm>

namespace MDR {

    // =========================================================================
    //  Optimized XORNegaBinaryBPEncoder
    //
    //  Key changes vs. the original:
    //
    //  1. last_2 / last_1 history is stored as flat arrays of T_stream words
    //     (one word per block per level), NOT vector<vector<vector<uchar>>>.
    //     This eliminates 3 levels of pointer chasing and keeps the data that
    //     decode_block touches in a single cache line per block.
    //
    //  2. The XOR recovery step operates on whole T_stream words with bitwise
    //     ops instead of looping over individual elements.
    //
    //  3. last_2 / last_1 update is also a whole-word copy — one assignment
    //     instead of a per-element loop.
    //
    //  4. The level table uses a flat vector instead of unordered_map.
    //
    //  5. Eliminated duplicated even/odd code paths by pre-computing a sign
    //     multiplier.
    // =========================================================================

    template<class T_data, class T_stream>
    class XORNegaBinaryBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        XORNegaBinaryBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "XORNegaBinaryBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "XORNegaBinaryBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "XORNegaBinaryBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "XORNegaBinaryBPEncoder: streams must be unsigned integers.");
        }

        // =====================================================================
        // encode (without level errors) — unchanged logic
        // =====================================================================
        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            exp += 2;
            const uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp  = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i = 0; i < num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            std::vector<T_stream *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }
            T_data const * data_pos = data;
            for(int32_t i = 0; i < n - (int32_t)block_size; i += block_size){
                for(uint32_t j = 0; j < block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = binary2negabinary((T_fps) shifted_data);
                }
                encode_block(int_data_buffer.data(), block_size, num_bitplanes, streams_pos);
            }
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j = 0; j < rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int_data_buffer[j] = binary2negabinary((T_fps) shifted_data);
                }
                encode_block(int_data_buffer.data(), rest_size, num_bitplanes, streams_pos);
            }
            for(int i = 0; i < num_bitplanes; i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            return streams;
        }

        // =====================================================================
        // encode (with level error collection) — unchanged logic
        // =====================================================================
        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            assert(num_bitplanes > 0);
            exp += 2;
            const uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp  = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i = 0; i < num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(n / UINT8_BITS + sizeof(T_stream)));
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            std::vector<T_stream *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }
            level_errors.clear();
            level_errors.resize(num_bitplanes + 1, 0.0);
            T_data const * data_pos = data;
            for(int32_t i = 0; i < n - (int32_t)block_size; i += block_size){
                for(uint32_t j = 0; j < block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_fps signed_int_data = (T_fps) shifted_data;
                    int_data_buffer[j] = binary2negabinary(signed_int_data);
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, shifted_data - signed_int_data, num_bitplanes);
                }
                encode_block(int_data_buffer.data(), block_size, num_bitplanes, streams_pos);
            }
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j = 0; j < rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_fps signed_int_data = (T_fps) shifted_data;
                    int_data_buffer[j] = binary2negabinary(signed_int_data);
                    collect_level_errors(level_errors, int_data_buffer[j], shifted_data, shifted_data - signed_int_data, num_bitplanes);
                }
                encode_block(int_data_buffer.data(), rest_size, num_bitplanes, streams_pos);
            }
            for(int i = 0; i < num_bitplanes; i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            for(size_t i = 0; i < level_errors.size(); i++){
                level_errors[i] = ldexp(level_errors[i], 2 * (-num_bitplanes + exp));
            }
            return streams;
        }

        // =====================================================================
        // decode — delegates to progressive_decode
        // =====================================================================
        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) {
            return progressive_decode(streams, n, exp, 0, num_bitplanes, streams.size());
        }

        // =====================================================================
        // progressive_decode — OPTIMIZED
        //
        //  • Flat contiguous history arrays (no vector<vector<vector>>)
        //  • Whole-word XOR recovery (no per-element loop for the XOR step)
        //  • Whole-word history shift (one assignment, not a loop)
        //  • Flat level table (vector, not unordered_map)
        //  • Merged even/odd code paths
        // =====================================================================
        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp,
                                     uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            const uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            const int32_t num_blocks = (n + block_size - 1) / block_size;

            using T_fps = typename std::conditional<std::is_same<T_data, double>::value, int64_t, int32_t>::type;
            using T_fp  = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;

            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            if(num_bitplanes == 0){
                memset(data, 0, n * sizeof(T_data));
                ensure_level_history(level, num_blocks);
                return data;
            }

            exp += 2;   // leave room for negabinary

            // --- Look up / create history for this level ---
            int decoded_bp = 0;
            ensure_level_history(level, num_blocks);
            if(level < (int)level_decoded_bp.size() && level_decoded_bp[level] > 0){
                decoded_bp = level_decoded_bp[level];
            }
            level_decoded_bp[level] += num_bitplanes;

            T_stream * hist2 = level_hist2[level].data();   // flat: num_blocks words
            T_stream * hist1 = level_hist1[level].data();

            // --- Set up stream cursors ---
            std::vector<T_stream const *> streams_pos(streams.size());
            for(size_t i = 0; i < streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream const *>(streams[i]);
            }

            // --- Decode ---
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            const T_data sign_mul = (ending_bitplane % 2 == 0) ? (T_data)1.0 : (T_data)-1.0;
            const T_data scale = ldexp(sign_mul, -(int)ending_bitplane + exp);

            std::vector<T_fp> int_data_buffer(block_size, 0);
            T_data * data_pos = data;

            const int remaining_first_bp = std::max(0, 2 - decoded_bp);

            for(int32_t bi = 0; bi < num_blocks; bi++){
                const int32_t base = bi * (int32_t)block_size;
                const int32_t cur_block_size = std::min((int32_t)block_size, n - base);

                memset(int_data_buffer.data(), 0, cur_block_size * sizeof(T_fp));

                T_stream h2 = hist2[bi];
                T_stream h1 = hist1[bi];

                // --- Decode block: bitplane loop with whole-word XOR ---
                for(int k = num_bitplanes - 1; k >= 0; k--){
                    const int bp_index = num_bitplanes - 1 - k;
                    T_stream bitplane_value = *(streams_pos[bp_index]++);

                    T_stream recovered;
                    if(k >= num_bitplanes - remaining_first_bp){
                        // First 1-2 bitplanes: no XOR, raw bits
                        recovered = bitplane_value;
                    } else {
                        // XOR recovery: whole-word operation
                        // Original per-element logic:
                        //   temp = last_2_bit ^ last_1_bit
                        //   recovered_bit = current_bit ? (1 - temp) : temp
                        // Equivalent whole-word:
                        //   temp = h2 ^ h1
                        //   recovered = bitplane_value ^ ~temp  (where current=1 → flip temp)
                        //             = bitplane_value ^ temp ^ all_ones... no.
                        //
                        // Let's derive carefully for the full block_size word:
                        //   temp = h2 ^ h1;
                        //   For bit i: recovered[i] = current[i] ? (1 - temp[i]) : temp[i]
                        //            = current[i] XOR temp[i]  ... wait:
                        //     current=0, temp=0 → 0  |  0 XOR 0 = 0 ✓
                        //     current=0, temp=1 → 1  |  0 XOR 1 = 1 ✓
                        //     current=1, temp=0 → 1  |  1 XOR 0 = 1 ✓
                        //     current=1, temp=1 → 0  |  1 XOR 1 = 0 ✓
                        // Yes! recovered = bitplane_value ^ (h2 ^ h1)
                        recovered = bitplane_value ^ (h2 ^ h1);
                    }

                    // For partial last block, mask off unused high bits
                    if(cur_block_size < (int32_t)block_size){
                        T_stream mask = (cur_block_size >= (int32_t)(sizeof(T_stream) * 8))
                                        ? ~(T_stream)0
                                        : ((T_stream)1 << cur_block_size) - 1;
                        recovered &= mask;
                    }

                    // Scatter recovered bits into int_data_buffer
                    for(int i = 0; i < cur_block_size; i++){
                        int_data_buffer[i] |= (T_fp)((recovered >> i) & 1u) << k;
                    }

                    // Update history: shift h2 ← h1, h1 ← recovered bitplane
                    h2 = h1;
                    h1 = recovered;
                }

                // Write back history for this block
                hist2[bi] = h2;
                hist1[bi] = h1;

                // Convert negabinary → binary → float
                for(int32_t j = 0; j < cur_block_size; j++){
                    *(data_pos++) = scale * (T_data) negabinary2binary(int_data_buffer[j]);
                }
            }

            return data;
        }

        void print() const {
            std::cout << "XORNegaBinary bitplane encoder" << std::endl;
        }

    private:
        // ----- helpers -----
        template<class T>
        uint32_t block_size_based_on_bitplane_int_type() const {
            if(std::is_same<T, uint64_t>::value) return 64;
            if(std::is_same<T, uint32_t>::value) return 32;
            if(std::is_same<T, uint16_t>::value) return 16;
            if(std::is_same<T, uint8_t>::value)  return 8;
            std::cerr << "Integer type not supported." << std::endl;
            exit(0);
        }

        inline uint64_t binary2negabinary(const int64_t x) const {
            return (x + (uint64_t)0xaaaaaaaaaaaaaaaaull) ^ (uint64_t)0xaaaaaaaaaaaaaaaaull;
        }
        inline uint32_t binary2negabinary(const int32_t x) const {
            return (x + (uint32_t)0xaaaaaaaau) ^ (uint32_t)0xaaaaaaaau;
        }
        inline int64_t negabinary2binary(const uint64_t x) const {
            return (x ^ 0xaaaaaaaaaaaaaaaaull) - 0xaaaaaaaaaaaaaaaaull;
        }
        inline int32_t negabinary2binary(const uint32_t x) const {
            return (x ^ 0xaaaaaaaau) - 0xaaaaaaaau;
        }

        inline void collect_level_errors(std::vector<double>& level_errors, uint32_t negabinary_data, float data, float mantissa, int num_bitplanes) const {
            level_errors[num_bitplanes] += mantissa * mantissa;
            for(int k = 1; k < num_bitplanes; k++){
                uint32_t mask = (1 << k) - 1;
                double diff = (double) negabinary2binary(negabinary_data & mask) + mantissa;
                level_errors[num_bitplanes - k] += diff * diff;
            }
            level_errors[0] += data * data;
        }

        template <class T_int>
        inline void encode_block(T_int const * data, size_t n, uint8_t num_bitplanes, std::vector<T_stream *>& streams_pos) const {
            for(int k = num_bitplanes - 1; k >= 0; k--){
                T_stream bitplane_value = 0;
                T_stream bitplane_index = num_bitplanes - 1 - k;
                if(k < num_bitplanes - 2){
                    for(size_t i = 0; i < n; i++){
                        T_stream last_2_bit = (T_stream)((data[i] >> (k + 2)) & 1u);
                        T_stream last_1_bit = (T_stream)((data[i] >> (k + 1)) & 1u);
                        T_stream current_bit = (T_stream)((data[i] >> k) & 1u);
                        bitplane_value += (T_stream)(last_2_bit ^ last_1_bit ^ current_bit) << i;
                    }
                } else {
                    for(size_t i = 0; i < n; i++){
                        bitplane_value += (T_stream)((data[i] >> k) & 1u) << i;
                    }
                }
                *(streams_pos[bitplane_index]++) = bitplane_value;
            }
        }

        // ----- Progressive decode state (flat arrays) -----

        // Per-level history: one T_stream word per block, storing the packed
        // bitplane values of the last-decoded and second-to-last-decoded
        // bitplanes.  This replaces vector<vector<vector<unsigned char>>>.
        mutable std::vector<std::vector<T_stream>> level_hist2;  // [level][block_idx]
        mutable std::vector<std::vector<T_stream>> level_hist1;
        mutable std::vector<int> level_decoded_bp;               // replaces unordered_map

        void ensure_level_history(int level, int32_t num_blocks) const {
            if(level >= (int)level_hist2.size()){
                level_hist2.resize(level + 1);
                level_hist1.resize(level + 1);
                level_decoded_bp.resize(level + 1, 0);
            }
            if((int)level_hist2[level].size() < num_blocks){
                level_hist2[level].assign(num_blocks, 0);
                level_hist1[level].assign(num_blocks, 0);
            }
        }
    };
}
#endif