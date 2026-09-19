#ifndef PRODM_ENCODER_BITPLANEENCODERINTERFACE_HPP
#define PRODM_ENCODER_BITPLANEENCODERINTERFACE_HPP

#include <vector>
#include <cstdint>
#include <cassert>

#include "ProDM/Namespace.hpp"

namespace ProDM {
    namespace concepts {
        #define UINT8_BITS 8 
        // concept of encoder which encodes T_data type data into bitstreams
        template<class T_data>
        class BitplaneEncoderInterface {
        public:

            virtual ~BitplaneEncoderInterface() = default;

            // virtual std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& streams_sizes) const = 0;

            virtual T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) = 0;

            virtual T_data * progressive_decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t starting_bitplane, uint8_t num_bitplanes, int level) = 0;

            virtual void print() const = 0;

        };
    }
}
#endif
