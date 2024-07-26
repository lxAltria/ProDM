#ifndef _MDR_WEIGHTED_PERBIT_BP_ENCODER_HPP
#define _MDR_WEIGHTED_PERBIT_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include "PerBitBPEncoder.hpp"
#include <bitset>
namespace MDR {

    #define PER_BIT_BLOCK_SIZE 1
    // per bit bitplane encoder that encodes data by bit using T_stream type buffer
    template<class T_data, class T_stream>
    class WeightedPerBitBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        WeightedPerBitBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            std::vector<uint8_t *> streams;
            exit(-1);
            return streams;
        }

        // normal encoding
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
            std::vector<BitEncoder> encoders;
            for(int i=0; i<streams.size(); i++){
                encoders.push_back(BitEncoder(reinterpret_cast<uint64_t*>(streams[i])));
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
                if(rest_size == 0) rest_size = block_size;
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

        // weighted encoding
        std::vector<uint8_t *> encode_weighted(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes, std::vector<double>& level_errors) const {
            // MGARD::writefile("refactor_coeff.dat", data, n);
            // if(n > 512*512*100) exit(0);
            // compute exp shift
            exp += num_weight_bitplanes;
            std::cout << "weighted encoding\n";
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
            std::vector<BitEncoder> encoders;
            for(int i=0; i<streams.size(); i++){
                encoders.push_back(BitEncoder(reinterpret_cast<uint64_t*>(streams[i])));
            }
            // init level errors
            level_errors.clear();
            level_errors.resize(num_bitplanes + 1);
            for(int i=0; i<level_errors.size(); i++){
                level_errors[i] = 0;
            }
            T_data const * data_pos = data;
            int * weights_pos = weights;
            for(int i=0; i<n - block_size; i+=block_size){
                T_stream sign_bitplane = 0;
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    // cur_data *= exp2(*(weights_pos++));
                    // T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp + *(weights_pos++));
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    // if(i+j==67989643){
                    //     std::cout << "exp = " << exp << ", num_bitplanes - exp = " << num_bitplanes - exp << std::endl;
                    //     std::cout << "data = " << data_pos[-1] << std::endl;
                    //     std::cout << "num_bitplanes = " << +num_bitplanes << std::endl;
                    //     std::cout << "weights = " << weights_pos[-1] << std::endl;
                    //     std::bitset<64> x(fp_data);
                    //     std::cout << x << std::endl;
                    //     for(int k=num_bitplanes - 1; k>=0; k--){
                    //         uint8_t bit = (fp_data >> k) & 1u;
                    //         std::cout << k << ": " << +bit << std::endl;
                    //     }
                    //     exit(0);                        
                    // }
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
                if(rest_size == 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    cur_data *= exp2(*(weights_pos++));
                    // T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp + *(weights_pos++));
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
            exit(-1);
            return NULL;
        }

        T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
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
            if(level_signs.size() == level){
                level_signs.push_back(std::vector<bool>(n, false));
                sign_flags.push_back(std::vector<bool>(n, false));
            }
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
                            uint8_t bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        // decode sign if possible
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            uint8_t bit = decoders[index].decode();
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
                            uint8_t bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            uint8_t bit = decoders[index].decode();
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

        T_data * progressive_decode_weighted(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) {
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            std::cout << "level_exp = " << exp << ", value = " << exp2(exp) << std::endl;
            exp += num_weight_bitplanes;
            std::cout << "weighted decoding\n";
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
            if(level_signs.size() == level){
                level_signs.push_back(std::vector<bool>(n, false));
                sign_flags.push_back(std::vector<bool>(n, false));
            }
            std::vector<bool>& signs = level_signs[level];
            std::vector<bool>& flags = sign_flags[level];
            const uint8_t ending_bitplane = starting_bitplane + num_bitplanes;
            // decode
            T_data * data_pos = data;
            int * weights_pos = weights;
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
                            uint8_t bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        // decode sign if possible
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            uint8_t bit = decoders[index].decode();
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
                    // T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
                    // cur_data /= exp2(*(weights_pos++));
                    T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp - *(weights_pos++));
                    *(data_pos++) = sign ? -cur_data : cur_data;
                    // if(i+j == 67989643){
                    //     std::cout << "starting_bitplane = " << +starting_bitplane << std::endl;
                    //     std::cout << "ending_bitplane = " << +ending_bitplane << std::endl;
                    //     std::cout << "exp = " << exp << ", - ending_bitplane + exp = " << - ending_bitplane + exp << std::endl;
                    //     std::cout << "data = " << data_pos[-1] << std::endl;
                    //     std::cout << "weights = " << weights_pos[-1] << std::endl;
                    //     std::bitset<64> x(fp_data);
                    //     std::cout << x << std::endl;
                    //     exit(0);                        
                    // }
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
                            uint8_t bit = decoders[index].decode();
                            fp_data += bit << k;
                        }
                    }
                    else{
                        bool first_bit = true;
                        for(int k=num_bitplanes - 1; k>=0; k--){
                            uint8_t index = num_bitplanes - 1 - k;
                            uint8_t bit = decoders[index].decode();
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
                    // T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp);
                    // cur_data /= exp2(*(weights_pos++));
                    T_data cur_data = ldexp((T_data)fp_data, - ending_bitplane + exp - *(weights_pos++));
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            // MGARD::writefile("reconstruct_coeff.dat", data, n);
            // if(n > 512*512*100) exit(0);
            return data;
        }
        void print() const {
            std::cout << "Per-bit bitplane encoder" << std::endl;
        }
        void set_weights(int * weights_, int num_weight_bitplanes_){
            weights = weights_;
            num_weight_bitplanes = num_weight_bitplanes_;
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
        int * weights;
        int num_weight_bitplanes = 0;
    };
}
#endif
