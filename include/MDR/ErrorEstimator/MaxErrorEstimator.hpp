#ifndef _MDR_MAX_ERROR_ESTIMATOR_HPP
#define _MDR_MAX_ERROR_ESTIMATOR_HPP

#include "ErrorEstimatorInterface.hpp"

namespace MDR {
    template<class T>
    class MaxErrorEstimator : public concepts::ErrorEstimatorInterface<T>{

    };
    // max error estimator for orthogonal basis
    template<class T>
    class MaxErrorEstimatorOB : public MaxErrorEstimator<T> {
    public:
        MaxErrorEstimatorOB(int num_dims){
            switch(num_dims){
                case 1:
                    c = 1.0 + sqrt(3)/2;
                    break;
                case 2:
                    c = 1.0 + 9.0/4;
                    break;
                case 3:
                    c = 1.0 + 21.0*sqrt(3)/8;
                    break;
                default:
                    std::cerr << num_dims << "-Dimentional error estimation not implemented." << std::endl;
                    exit(-1);
            }
	    c *= 4; // 2 more bitplane for negabinary
        }
        MaxErrorEstimatorOB() : MaxErrorEstimatorOB(1) {}

        inline T estimate_error(T error, int level) const {
            return c * error;
        }
        /////////
        inline T estimate_error(T error, int level, int num_levels) const {
            return c * error;
        }
        inline T estimate_error(T data, T reconstructed_data, int level) const {
            return c * (data - reconstructed_data);
        }
        inline T estimate_error_gain(T base, T current_level_err, T next_level_err, int level) const {
            return c * (current_level_err - next_level_err);
        }
        void print() const {
            std::cout << "Max absolute error estimator (up to 3 dimensions) for orthogonal basis." << std::endl;
        }
    private:
        // derived constant
        T c = 0;
    };
    // max error estimator for hierarchical basis
    // c = 1 as all the operations are linear
    template<class T>
    class MaxErrorEstimatorHB : public MaxErrorEstimator<T> {
    public:
        MaxErrorEstimatorHB(){}
        inline T estimate_error(T error, int level) const {
            return error;
        }
        ///
        inline T estimate_error(T error, int level, int num_levels) const {
            return error;
        }
        inline T estimate_error(T data, T reconstructed_data, int level) const {
            return data - reconstructed_data;
        }
        inline T estimate_error_gain(T base, T current_level_err, T next_level_err, int level) const {
            return current_level_err - next_level_err;
        }
        void print() const {
            std::cout << "Max absolute error estimator for hierarchical basis." << std::endl;
        }
    };
    template<class T>
    class MaxErrorEstimatorHB_new : public MaxErrorEstimator<T> {
    public:
        MaxErrorEstimatorHB_new(int num_dims){
            switch(num_dims){
                case 1:
                    c = 1;
                    break; 
                case 2:
                    c = 2;
                    break;
                case 3:
                    c = 3;
                    break;
                default:
                    std::cerr << num_dims << "-Dimentional error estimation not implemented." << std::endl;
                    exit(-1);
            }
        }
        MaxErrorEstimatorHB_new() : MaxErrorEstimatorHB_new(1) {}
        inline T estimate_error(T error, int level) const {
            return (level) ? c * error : error;
        }
        ///
        inline T estimate_error(T error, int level, int num_levels) const {
            return (level) ? c * error : error;;
        }
        inline T estimate_error(T data, T reconstructed_data, int level) const {
            return (level) ? c * (data - reconstructed_data) : data - reconstructed_data;
        }
        inline T estimate_error_gain(T base, T current_level_err, T next_level_err, int level) const {
            return (level) ? c * (current_level_err - next_level_err) : current_level_err - next_level_err;
        }
        void print() const {
            std::cout << "Max absolute error estimator for hierarchical basis." << std::endl;
        }
    private:
        T c = 0;
    };
    template<class T>
    class MaxErrorEstimatorHBCubic : public MaxErrorEstimator<T> {
    public:
        MaxErrorEstimatorHBCubic(int num_dims){
            switch(num_dims){
                case 1:
                    c = 5.0/4;
                    break;
                case 2:
                    c = (5.0/4) * (5.0/4);
                    break;
                case 3:
                    c = (5.0/4) * (5.0/4) * (5.0/4);
                    break;
                default:
                    std::cerr << num_dims << "-Dimentional error estimation not implemented." << std::endl;
                    exit(-1);
            }
        }
        MaxErrorEstimatorHBCubic() : MaxErrorEstimatorHBCubic(1) {}
        inline T estimate_error(T error, int level) const {
            return c*error;
        }
        inline T estimate_error(T error, int level, int num_levels) const {
            if(first_time){
                C = std::vector<T>(num_levels, 1);
                for(int i=num_levels-2; i>=0; i--){
                    C[i] = C[i+1] * c;
                }
                first_time = false;
            }
            return C[level] * error;
        }
        inline T estimate_error(T data, T reconstructed_data, int level) const {
            return c * (data - reconstructed_data);
        }
        inline T estimate_error_gain(T base, T current_level_err, T next_level_err, int level) const {
            return c * (current_level_err - next_level_err);
        }
        void print() const {
            std::cout << "Max absolute error estimator for hierarchical basis." << std::endl;
        }
    private:
        double c = 0;
        mutable std::vector<T> C;
        mutable bool first_time = true;
    };    
}
#endif
