#ifndef PRODM_ERRORCONTROL_SIZEINTERPRETER_SIZEINTERPRETERINTERFACE_HPP
#define PRODM_ERRORCONTROL_SIZEINTERPRETER_SIZEINTERPRETERINTERFACE_HPP

#include <vector>
#include <cstdint>

#include "ProDM/Namespace.hpp"

namespace ProDM {
    namespace concepts {

        // level bit-plane reorganizer: EBCOT-like algorithm for multilevel bit-plane truncation
        class SizeInterpreterInterface {
        public:

            virtual ~SizeInterpreterInterface() = default;

            virtual std::vector<uint32_t> interpret_retrieve_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
