#ifndef _MDR_TUNER_INTERFACE_HPP
#define _MDR_TUNER_INTERFACE_HPP
#include <cstdint>

namespace MDR{
    namespace concepts{

        // Tuner
        template<class T>
        class TunerInterface {
        public:

            virtual ~TunerInterface() = default;

            virtual void tune(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes, uint32_t stride, uint32_t block_size) = 0;

            virtual void print() const = 0;
        };
    }
}

#endif