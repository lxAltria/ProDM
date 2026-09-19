#ifndef PRODM_ERRORCONTROL_COLLECTOR_ERRORCOLLECTORINTERFACE_HPP
#define PRODM_ERRORCONTROL_COLLECTOR_ERRORCOLLECTORINTERFACE_HPP

#include <vector>
#include <cstdint>

#include "ProDM/Namespace.hpp"

namespace ProDM {
    namespace concepts {

        // Error estimator: estimate impact of level errors on the final error
        template<class T>
        class ErrorCollectorInterface {
        public:

            virtual ~ErrorCollectorInterface() = default;

            virtual std::vector<double> collect_level_error(T const * data, size_t n, int num_bitplanes, T max_level_error) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
