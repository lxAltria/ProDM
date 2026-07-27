#ifndef _PDR_RECONSTRUCTOR_INTERFACE_HPP
#define _PDR_RECONSTRUCTOR_INTERFACE_HPP
#include <cstdint>

namespace PDR {
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
