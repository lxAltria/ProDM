#ifndef _PDR_APPROXIMATOR_INTERFACE_HPP
#define _PDR_APPROXIMATOR_INTERFACE_HPP

namespace PDR {
    namespace concepts {

        // inplace data approximator: de-correlates and overwrites original data
        template<class T>
        class ApproximatorInterface {
        public:

            virtual ~ApproximatorInterface() = default;

            virtual void refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximation_eb, std::vector<uint32_t> strides) = 0;

            virtual void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::vector<uint32_t> strides) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
