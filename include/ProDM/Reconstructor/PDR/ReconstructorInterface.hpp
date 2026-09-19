#ifndef PRODM_RECONSTRUCTOR_PDR_RECONSTRUCTORINTERFACE_HPP
#define PRODM_RECONSTRUCTOR_PDR_RECONSTRUCTORINTERFACE_HPP
#include <cstdint>
#include <vector>

#include "ProDM/Namespace.hpp"

namespace ProDM::PDR {
    namespace concepts {

        // reconstructor: a general interface for scientific data reconstructor
        template<class T>
        class ReconstructorInterface {
        public:

            virtual ~ReconstructorInterface() = default;

            virtual T * reconstruct(double tolerance) = 0;

            virtual T * progressive_reconstruct(double tolerance) = 0;

            virtual void load_metadata() = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
