#ifndef _PDR_DUMMY_APPROXIMATOR_HPP
#define _PDR_DUMMY_APPROXIMATOR_HPP

#include "ApproximatorInterface.hpp"

namespace PDR {
    // Dummy approximator with 0 prediction
    template<class T>
    class DummyApproximator : public concepts::ApproximatorInterface<T> {
    public:
        DummyApproximator(){}
        size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximator_eb, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            return 0;
        }

        void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
        	return;
        }

        size_t get_size() const {
            return 0;
        }

        void print() const {
            std::cout << "PDR dummy approximator" << std::endl;
        }
    };
}

#endif