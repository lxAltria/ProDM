#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Legacy/QoIUtils.hpp"
#include "ProDM/Utils/MaskUtils.hpp"
#include "ProDM/Reconstructor/PDR/Reconstructor.hpp"
#include "ProDM/Utils/StatUtils.hpp"
#define BP 0
#define WBP 1
#define QoI_Vtot 0
#define QoI_Vtot2 1

using namespace std;
using namespace ProDM;
using namespace ProDM::Legacy;

const vector<string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};

// Error-control iterations: estimate the QoI error from the per-variable error
// bounds and decrease the bounds until the QoI tolerance is met.
// Each QoI has its own set of functions since QoIs are computed in different ways;
// the V_TOT functions follow test/qoi_Vtot_d64.cpp and the V_TOT_2 functions
// follow test/qoi_Vtot2_d64.cpp of the original evaluation code.

// ------------------------------- QoI: V_TOT --------------------------------
template<class T>
bool halving_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const T * V_TOT_ori, std::vector<T>& error_V_TOT, std::vector<T>& error_est_V_TOT){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        // error of total velocity square
        T e_V_TOT_2 = 0;
        if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        // error of total velocity
        T e_V_TOT = 0;
        if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
        T V_TOT = sqrt(V_TOT_2);

        error_est_V_TOT[i] = e_V_TOT;
        error_V_TOT[i] = V_TOT - V_TOT_ori[i];

        if(max_value < error_est_V_TOT[i]){
            max_value = error_est_V_TOT[i];
            max_index = i;
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            eb_Vx = eb_Vx / 1.5;
            eb_Vy = eb_Vy / 1.5;
            eb_Vz = eb_Vz / 1.5;
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
            estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const T * V_TOT_ori, std::vector<T>& error_V_TOT, std::vector<T>& error_est_V_TOT){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        // error of total velocity square
        T e_V_TOT_2 = 0;
        if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        // error of total velocity
        T e_V_TOT = 0;
        if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
        T V_TOT = sqrt(V_TOT_2);

        error_est_V_TOT[i] = e_V_TOT;
        error_V_TOT[i] = V_TOT - V_TOT_ori[i];

        if(max_value < error_est_V_TOT[i]){
            max_value = error_est_V_TOT[i];
            max_index = i;
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            T estimate_error_Vx = 0;
            {
                T eb_Vx_ = eb_Vx / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx_) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
                estimate_error_Vx = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            T estimate_error_Vy = 0;
            {
                T eb_Vy_ = eb_Vy / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy_) + compute_bound_x_square(Vz[i], eb_Vz);
                estimate_error_Vy = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            T estimate_error_Vz = 0;
            {
                T eb_Vz_ = eb_Vz / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz_);
                estimate_error_Vz = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            const T epsilon = 1e-6;
            T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz});
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
            if (ebs[0] / eb_Vx > 10 || ebs[1] / eb_Vy > 10 || ebs[2] / eb_Vz > 10) break;
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const std::vector<std::vector<int>>& weights, const T * V_TOT_ori, std::vector<T>& error_V_TOT, std::vector<T>& error_est_V_TOT){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        if(mask[i]) {
            // error of total velocity square
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
            // error of total velocity
            T e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            T V_TOT = sqrt(V_TOT_2);

            error_est_V_TOT[i] = e_V_TOT;
            error_V_TOT[i] = V_TOT - V_TOT_ori[i];

            if(max_value < error_est_V_TOT[i]){
                max_value = error_est_V_TOT[i];
                max_index = i;
            }
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            eb_Vx = eb_Vx / 1.5;
            eb_Vy = eb_Vy / 1.5;
            eb_Vz = eb_Vz / 1.5;
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const std::vector<std::vector<int>>& weights, const T * V_TOT_ori, std::vector<T>& error_V_TOT, std::vector<T>& error_est_V_TOT){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        if(mask[i]) {
            // error of total velocity square
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
            // error of total velocity
            T e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            T V_TOT = sqrt(V_TOT_2);

            error_est_V_TOT[i] = e_V_TOT;
            error_V_TOT[i] = V_TOT - V_TOT_ori[i];

            if(max_value < error_est_V_TOT[i]){
                max_value = error_est_V_TOT[i];
                max_index = i;
            }
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            T estimate_error_Vx = 0;
            {
                T eb_Vx_ = eb_Vx / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx_, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
                estimate_error_Vx = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            T estimate_error_Vy = 0;
            {
                T eb_Vy_ = eb_Vy / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy_, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
                estimate_error_Vy = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            T estimate_error_Vz = 0;
            {
                T eb_Vz_ = eb_Vz / 1.5;
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz_, -weights[2][i]));
                estimate_error_Vz = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
            T sum_err = 3 * estimate_error - (estimate_error_Vx + estimate_error_Vy + estimate_error_Vz);
            T w_Vx = (estimate_error - estimate_error_Vx) / sum_err;
            T w_Vy = (estimate_error - estimate_error_Vy) / sum_err;
            T w_Vz = (estimate_error - estimate_error_Vz) / sum_err;
            // smooth proportional update
            T factor_base = 1.5;
            T alpha = 1.0 - (1.0 / factor_base);
            eb_Vx *= (1.0 - alpha * w_Vx);
            eb_Vy *= (1.0 - alpha * w_Vy);
            eb_Vz *= (1.0 - alpha * w_Vz);
            {
                T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
                estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
            }
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

// ------------------------------ QoI: V_TOT_2 -------------------------------
template<class T>
bool halving_error_V_TOT2_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const T * V_TOT2_ori, std::vector<T>& error_V_TOT2, std::vector<T>& error_est_V_TOT2){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        // error of total velocity square
        T e_V_TOT_2 = 0;
        if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];

        error_est_V_TOT2[i] = e_V_TOT_2;
        error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

        if(max_value < error_est_V_TOT2[i]){
            max_value = error_est_V_TOT2[i];
            max_index = i;
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            eb_Vx = eb_Vx / 1.5;
            eb_Vy = eb_Vy / 1.5;
            eb_Vz = eb_Vz / 1.5;
            estimate_error = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT2_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const T * V_TOT2_ori, std::vector<T>& error_V_TOT2, std::vector<T>& error_est_V_TOT2){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        // error of total velocity square
        T e_V_TOT_2 = 0;
        if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
        T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];

        error_est_V_TOT2[i] = e_V_TOT_2;
        error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

        if(max_value < error_est_V_TOT2[i]){
            max_value = error_est_V_TOT2[i];
            max_index = i;
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            T estimate_error_Vx = 0;
            {
                T eb_Vx_ = eb_Vx / 1.5;
                estimate_error_Vx = compute_bound_x_square(Vx[i], eb_Vx_) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
            }
            T estimate_error_Vy = 0;
            {
                T eb_Vy_ = eb_Vy / 1.5;
                estimate_error_Vy = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy_) + compute_bound_x_square(Vz[i], eb_Vz);
            }
            T estimate_error_Vz = 0;
            {
                T eb_Vz_ = eb_Vz / 1.5;
                estimate_error_Vz = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz_);
            }
            const T epsilon = 1e-6;
            T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz});
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
            if (ebs[0] / eb_Vx > 10 || ebs[1] / eb_Vy > 10 || ebs[2] / eb_Vz > 10) break;
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT2_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const std::vector<std::vector<int>>& weights, const T * V_TOT2_ori, std::vector<T>& error_V_TOT2, std::vector<T>& error_est_V_TOT2){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        if(mask[i]) {
            // error of total velocity square
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];

            error_est_V_TOT2[i] = e_V_TOT_2;
            error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

            if(max_value < error_est_V_TOT2[i]){
                max_value = error_est_V_TOT2[i];
                max_index = i;
            }
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            eb_Vx = eb_Vx / 1.5;
            eb_Vy = eb_Vy / 1.5;
            eb_Vz = eb_Vz / 1.5;
            estimate_error = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template<class T>
bool halving_error_V_TOT2_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, const std::vector<std::vector<int>>& weights, const T * V_TOT2_ori, std::vector<T>& error_V_TOT2, std::vector<T>& error_est_V_TOT2){
    T eb_Vx = ebs[0];
    T eb_Vy = ebs[1];
    T eb_Vz = ebs[2];
    T max_value = 0;
    size_t max_index = 0;
    for(size_t i=0; i<n; i++){
        if(mask[i]) {
            // error of total velocity square
            T e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];

            error_est_V_TOT2[i] = e_V_TOT_2;
            error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

            if(max_value < error_est_V_TOT2[i]){
                max_value = error_est_V_TOT2[i];
                max_index = i;
            }
        }
    }
    if(max_value > tau){
        auto i = max_index;
        T estimate_error = max_value;
        T eb_Vx = ebs[0];
        T eb_Vy = ebs[1];
        T eb_Vz = ebs[2];
        while(estimate_error > tau){
            T estimate_error_Vx = 0;
            {
                T eb_Vx_ = eb_Vx / 1.5;
                estimate_error_Vx = compute_bound_x_square(Vx[i], ldexp(eb_Vx_, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            }
            T estimate_error_Vy = 0;
            {
                T eb_Vy_ = eb_Vy / 1.5;
                estimate_error_Vy = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy_, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            }
            T estimate_error_Vz = 0;
            {
                T eb_Vz_ = eb_Vz / 1.5;
                estimate_error_Vz = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz_, -weights[2][i]));
            }
            T sum_err = 3 * estimate_error - (estimate_error_Vx + estimate_error_Vy + estimate_error_Vz);
            T w_Vx = (estimate_error - estimate_error_Vx) / sum_err;
            T w_Vy = (estimate_error - estimate_error_Vy) / sum_err;
            T w_Vz = (estimate_error - estimate_error_Vz) / sum_err;
            // smooth proportional update
            T factor_base = 1.5;
            T alpha = 1.0 - (1.0 / factor_base);
            eb_Vx *= (1.0 - alpha * w_Vx);
            eb_Vy *= (1.0 - alpha * w_Vy);
            eb_Vz *= (1.0 - alpha * w_Vz);
            {
                estimate_error = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
            }
        }
        ebs[0] = eb_Vx;
        ebs[1] = eb_Vy;
        ebs[2] = eb_Vz;
        return false;
    }
    return true;
}

template <class T>
void launch_reconstructor(const string& data_dir, const string& refactor_dir, const vector<double>& tolerance, const string& suffix, int mode, int qoi, int decrease_method){
#ifdef PRODM_HAVE_HPEZ
    using T_stream = uint32_t;
    int n_variable = var_list.size();
    size_t num_elements = 0;
    vector<vector<T>> vars_ori(n_variable);
    for(int i=0; i<n_variable; i++){
        string filename = data_dir + "/" + var_list[i] + suffix;
        vars_ori[i] = ProDM::readfile<T>(filename.c_str(), num_elements);
        if(num_elements == 0){
            cerr << "Cannot read " << filename << endl;
            exit(-1);
        }
    }

    // base error bound: relative to the value range across the three velocities
    vector<T> V3(3 * num_elements);
    for(int i=0; i<n_variable; i++){
        memcpy(V3.data() + i * num_elements, vars_ori[i].data(), num_elements * sizeof(T));
    }
    T eb_base = compute_value_range(V3);

    // original QoI values; the tolerances are relative to the QoI value range
    string qoi_name = (qoi == QoI_Vtot) ? "V_TOT" : "V_TOT_2";
    vector<T> QoI_ori(num_elements);
    if(qoi == QoI_Vtot) compute_VTOT(vars_ori[0].data(), vars_ori[1].data(), vars_ori[2].data(), num_elements, QoI_ori.data());
    else compute_VTOT2(vars_ori[0].data(), vars_ori[1].data(), vars_ori[2].data(), num_elements, QoI_ori.data());
    T tau_base = compute_value_range(QoI_ori);

    string mask_file = refactor_dir + "/mask.bin";
    uint32_t mask_file_size = 0;
    auto mask = ProDM::readmask(mask_file.c_str(), mask_file_size);
    if(mask.size() != num_elements){
        cerr << "Mask " << mask_file << " has " << mask.size() << " elements, expected " << num_elements << "; run the refactor first" << endl;
        exit(-1);
    }

    const int max_iter = 30;
    vector<vector<T>> reconstructed_vars(n_variable, vector<T>(num_elements));
    vector<T> error_QoI(num_elements);
    vector<T> error_est_QoI(num_elements);

    if(mode == BP){
        vector<PDR::ApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, ProDM::NegaBinaryBPEncoder<T, T_stream>, ProDM::AdaptiveLevelCompressor, ProDM::SignExcludeGreedyBasedSizeInterpreter<ProDM::MaxErrorEstimatorHB<T>>, ProDM::MaxErrorEstimatorHB<T>, ProDM::ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            string rdir_prefix = refactor_dir + "/" + var_list[i] + "_refactored";
            string metadata_file = rdir_prefix + "/metadata.bin";
            vector<string> files = {rdir_prefix + "/level_0.bin"};
            auto approximator = PDR::HPEZApproximator<T>();
            auto encoder = ProDM::NegaBinaryBPEncoder<T, T_stream>();
            auto compressor = ProDM::AdaptiveLevelCompressor(64);
            auto estimator = ProDM::MaxErrorEstimatorHB<T>();
            auto interpreter = ProDM::SignExcludeGreedyBasedSizeInterpreter<ProDM::MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(PDR::ApproximationBasedReconstructor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>(approximator, encoder, compressor, interpreter, retriever));
            reconstructors.back().mask = mask;
            reconstructors.back().load_metadata();
        }
        for(int t=0; t<tolerance.size(); t++){
            T tau = tau_base * tolerance[t];
            vector<T> ebs(n_variable, eb_base * tolerance[t]);
            int iter = 0;
            bool tolerance_met = false;
            T max_act_error = 0, max_est_error = 0;
            vector<size_t> total_retrieved_size(n_variable, 0);
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            while((!tolerance_met) && (iter < max_iter)){
                iter++;
                for(int i=0; i<n_variable; i++){
                    auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                    total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
                if(qoi == QoI_Vtot){
                    if(!decrease_method) tolerance_met = halving_error_V_TOT_uniform(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, QoI_ori.data(), error_QoI, error_est_QoI);
                    else tolerance_met = halving_error_V_TOT_coordinate(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, QoI_ori.data(), error_QoI, error_est_QoI);
                }
                else{
                    if(!decrease_method) tolerance_met = halving_error_V_TOT2_uniform(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, QoI_ori.data(), error_QoI, error_est_QoI);
                    else tolerance_met = halving_error_V_TOT2_coordinate(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, QoI_ori.data(), error_QoI, error_est_QoI);
                }
                max_act_error = print_max_abs(qoi_name + " error", error_QoI);
                max_est_error = print_max_abs(qoi_name + " error_est", error_est_QoI);
            }
            clock_gettime(CLOCK_REALTIME, &end);
            double elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
            size_t total_size = mask_file_size + std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), size_t(0));
            cout << "requested_error = " << tau << endl;
            cout << "max_est_error = " << max_est_error << endl;
            cout << "max_act_error = " << max_act_error << endl;
            cout << "iter = " << iter << endl;
            cout << "aggregated cr = " << n_variable * num_elements * sizeof(T) * 1.0 / total_size << endl;
            cout << "bitrate = " << total_size * 8.0 / (n_variable * num_elements) << endl;
            printf("elapsed_time = %.6f\n", elapsed_time);
            cout << endl;
        }
    }
    else{
        vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, ProDM::WeightedNegaBinaryBPEncoder<T, T_stream>, ProDM::AdaptiveLevelCompressor, ProDM::SignExcludeGreedyBasedSizeInterpreter<ProDM::MaxErrorEstimatorHB<T>>, ProDM::MaxErrorEstimatorHB<T>, ProDM::ConcatLevelFileRetriever>> reconstructors;
        vector<vector<int>> weights(n_variable, vector<int>(num_elements));
        for(int i=0; i<n_variable; i++){
            string rdir_prefix = refactor_dir + "/" + var_list[i] + "_refactored";
            string metadata_file = rdir_prefix + "/metadata.bin";
            vector<string> files = {rdir_prefix + "/level_0.bin"};
            auto approximator = PDR::HPEZApproximator<T>();
            auto encoder = ProDM::WeightedNegaBinaryBPEncoder<T, T_stream>();
            auto compressor = ProDM::AdaptiveLevelCompressor(64);
            auto estimator = ProDM::MaxErrorEstimatorHB<T>();
            auto interpreter = ProDM::SignExcludeGreedyBasedSizeInterpreter<ProDM::MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(PDR::WeightedApproximationBasedReconstructor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>(approximator, encoder, compressor, interpreter, retriever));
            reconstructors.back().mask = mask;
            // the weights are stored with the first variable and shared by the others
            if(i == 0) reconstructors.back().fetch_weight = true;
            else reconstructors.back().copy_int_weights(weights[0], reconstructors[0].get_max_weight());
            reconstructors.back().load_metadata();
            weights[i] = reconstructors.back().get_int_weights();
        }
        size_t weight_file_size = reconstructors[0].get_weight_file_size();
        for(int t=0; t<tolerance.size(); t++){
            T tau = tau_base * tolerance[t];
            vector<T> ebs(n_variable);
            for(int i=0; i<n_variable; i++){
                ebs[i] = ldexp(eb_base * (T)tolerance[t], reconstructors[0].get_max_weight());
            }
            int iter = 0;
            bool tolerance_met = false;
            T max_act_error = 0, max_est_error = 0;
            vector<size_t> total_retrieved_size(n_variable, 0);
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            while((!tolerance_met) && (iter < max_iter)){
                iter++;
                for(int i=0; i<n_variable; i++){
                    auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
                    total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
                if(qoi == QoI_Vtot){
                    if(!decrease_method) tolerance_met = halving_error_V_TOT_uniform(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, weights, QoI_ori.data(), error_QoI, error_est_QoI);
                    else tolerance_met = halving_error_V_TOT_coordinate(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, weights, QoI_ori.data(), error_QoI, error_est_QoI);
                }
                else{
                    if(!decrease_method) tolerance_met = halving_error_V_TOT2_uniform(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, weights, QoI_ori.data(), error_QoI, error_est_QoI);
                    else tolerance_met = halving_error_V_TOT2_coordinate(reconstructed_vars[0].data(), reconstructed_vars[1].data(), reconstructed_vars[2].data(), num_elements, mask, tau, ebs, weights, QoI_ori.data(), error_QoI, error_est_QoI);
                }
                max_act_error = print_max_abs(qoi_name + " error", error_QoI);
                max_est_error = print_max_abs(qoi_name + " error_est", error_est_QoI);
            }
            clock_gettime(CLOCK_REALTIME, &end);
            double elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
            size_t total_size = mask_file_size + weight_file_size + std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), size_t(0));
            cout << "requested_error = " << tau << endl;
            cout << "max_est_error = " << max_est_error << endl;
            cout << "max_act_error = " << max_act_error << endl;
            cout << "iter = " << iter << endl;
            cout << "aggregated cr = " << n_variable * num_elements * sizeof(T) * 1.0 / total_size << endl;
            cout << "bitrate = " << total_size * 8.0 / (n_variable * num_elements) << endl;
            printf("elapsed_time = %.6f\n", elapsed_time);
            cout << endl;
        }
    }
#else
    cerr << "This tool requires the HPEZ approximator; rebuild with -DPRODM_WITH_HPEZ=ON" << endl;
    exit(-1);
#endif
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_dir refactor_dir num_tolerance tolerance1 ... toleranceN -[dataType: f/d] [mode: BP-0, WBP-1] [QoI: Vtot-0, Vtot2-1] [decrease_method: uniform-0, coordinate-1]"
                  << std::endl
                  << "tolerances are relative error bounds on the QoI; reads VelocityX/Y/Z (.dat for -d, .dat.f32 for -f) from data_dir" << std::endl
                  << "example: " << cmd <<
                  " data refactor 3 1e-1 1e-2 1e-3 -d 1 0 0" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    std::string data_dir = string(argv[argv_id++]);
    std::string refactor_dir = string(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id++]);
    if(num_tolerance <= 0 || argc < argv_id + num_tolerance + 3){
        std::cerr << "Insufficient or invalid arguments (num_tolerance parsed as " << num_tolerance << "); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id++]);
    }

    std::string dtype = string(argv[argv_id++]);
    int mode = atoi(argv[argv_id++]);
    int qoi = atoi(argv[argv_id++]);
    if((mode != BP && mode != WBP) || (qoi != QoI_Vtot && qoi != QoI_Vtot2)){
        std::cerr << "Unknown mode " << mode << " (expected BP-0 or WBP-1) or QoI " << qoi << " (expected Vtot-0 or Vtot2-1)" << std::endl;
        usage(argv[0]);
        return -1;
    }
    int decrease_method = 0;
    if(argc > argv_id) decrease_method = atoi(argv[argv_id++]);
    if(decrease_method != 0 && decrease_method != 1){
        std::cerr << "Unknown decrease_method " << decrease_method << " (expected uniform-0 or coordinate-1)" << std::endl;
        usage(argv[0]);
        return -1;
    }

    if (strcmp(dtype.c_str(), "-f") == 0){
        launch_reconstructor<float>(data_dir, refactor_dir, tolerance, ".dat.f32", mode, qoi, decrease_method);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        launch_reconstructor<double>(data_dir, refactor_dir, tolerance, ".dat", mode, qoi, decrease_method);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }

    return 0;
}
