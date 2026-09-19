#ifndef PRODM_DECOMPOSER_APPROXIMATION_APPROXIMATORINTERFACE_HPP
#define PRODM_DECOMPOSER_APPROXIMATION_APPROXIMATORINTERFACE_HPP

#include <string>

#include <vector>
#include <cstdint>
#include <cstddef>

#include "ProDM/Namespace.hpp"

namespace ProDM::PDR {
    namespace concepts {

        // inplace data approximator: de-correlates and overwrites original data
        template<class T>
        class ApproximatorInterface {
        public:

            virtual ~ApproximatorInterface() = default;

            virtual size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximation_eb, std::string filename, std::vector<uint32_t> strides) = 0;

            virtual void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename, std::vector<uint32_t> strides) = 0;

            virtual size_t get_size() const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
