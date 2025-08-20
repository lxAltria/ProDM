#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include<limits>
#include "mpi.h"
#include <sstream>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "nomask_Synthesizer4GE.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#define Dummy 0
#define SZ3 1
#define PMGARD 2
#define GE 3
using namespace MDR;

using T = double;
std::vector<double> P_ori;
std::vector<double> D_ori;
std::vector<double> Vx_ori;
std::vector<double> Vy_ori;
std::vector<double> Vz_ori;
double * P_dec = NULL;
double * D_dec = NULL;
double * Vx_dec = NULL;
double * Vy_dec = NULL;
double * Vz_dec = NULL;
double * Mach_ori = NULL;
std::vector<double> error_Mach;
std::vector<double> error_est_Mach;
int iter = 0;
double local_elapsed_time = 0;

template<class T>
bool halfing_error_Mach_uniform(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T eb_P = ebs[3];
	T eb_D = ebs[4];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	T max_value = 0;
	int max_index = 0;
    T max_Vx = 0;
    T max_Vy = 0;
    T max_Vz = 0;
    T max_P = 0;
    T max_D = 0;
    T max_e_VTOT_2 = 0;
    T max_VTOT_2 = 0;
    T max_e_T = 0;
    T max_T = 0;
    T max_C = 0;
    T max_e_C = 0;
    T max_e_Mach = 0;
    T max_Mach = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);
		// error of C
		T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		T C = c_2 * sqrt(Temp);
		// print_error("C", C, C_ori[i], e_C);
		T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		T Mach = V_TOT / C;
		// print_error("Mach", Mach, Mach_ori[i], e_Mach);
		error_est_Mach[i] = e_Mach;
		error_Mach[i] = Mach - Mach_ori[i];
		if(max_value < error_est_Mach[i]){
			max_value = error_est_Mach[i];
			max_index = i;
            max_Vx = Vx[i];
            max_Vy = Vy[i];
            max_Vz = Vz[i];
            max_P = P[i];
            max_D = D[i];
            max_e_VTOT_2 = e_V_TOT_2;
            max_VTOT_2 = V_TOT_2;
            max_e_T = e_T;
            max_T = Temp;
            max_e_C = e_C;
            max_C = C;
            max_e_Mach = e_Mach;
            max_Mach = Mach;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
    // std::cout << "e_Mach = " << max_e_Mach << ", Mach = " << max_Mach << std::endl;
    // std::cout << "e_C = " << max_e_C << ", C = " << max_C << std::endl;
    // std::cout << "e_VTOT_2 = " << max_e_VTOT_2 << ", VTOT_2 = " << max_VTOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
    // std::cout << "e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		T eb_P = ebs[3];
		T eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			T e_V_TOT_2 = 0;
			if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
			T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			T e_V_TOT = 0;
			if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			T V_TOT = sqrt(V_TOT_2);
			T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			T Temp = P[i] / (D[i] * R);
			T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
			T C = c_2 * sqrt(Temp);
			estimate_error = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
            if (ebs[0] / eb_Vx > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}
    return true;
}

template<class T>
bool halfing_error_Mach_coordinate(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T eb_P = ebs[3];
	T eb_D = ebs[4];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	T max_value = 0;
	int max_index = 0;
    T max_Vx = 0;
    T max_Vy = 0;
    T max_Vz = 0;
    T max_P = 0;
    T max_D = 0;
    T max_e_VTOT_2 = 0;
    T max_VTOT_2 = 0;
    T max_e_T = 0;
    T max_T = 0;
    T max_C = 0;
    T max_e_C = 0;
    T max_e_Mach = 0;
    T max_Mach = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);
		// error of C
		T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		T C = c_2 * sqrt(Temp);
		// print_error("C", C, C_ori[i], e_C);
		T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		T Mach = V_TOT / C;
		// print_error("Mach", Mach, Mach_ori[i], e_Mach);
		error_est_Mach[i] = e_Mach;
		error_Mach[i] = Mach - Mach_ori[i];
		if(max_value < error_est_Mach[i]){
			max_value = error_est_Mach[i];
			max_index = i;
            max_Vx = Vx[i];
            max_Vy = Vy[i];
            max_Vz = Vz[i];
            max_P = P[i];
            max_D = D[i];
            max_e_VTOT_2 = e_V_TOT_2;
            max_VTOT_2 = V_TOT_2;
            max_e_T = e_T;
            max_T = Temp;
            max_e_C = e_C;
            max_C = C;
            max_e_Mach = e_Mach;
            max_Mach = Mach;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
    // std::cout << "e_Mach = " << max_e_Mach << ", Mach = " << max_Mach << std::endl;
    // std::cout << "e_C = " << max_e_C << ", C = " << max_C << std::endl;
    // std::cout << "e_VTOT_2 = " << max_e_VTOT_2 << ", VTOT_2 = " << max_VTOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
    // std::cout << "e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    
    // ************************************************ Vx, Vy, Vz, P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		T eb_P = ebs[3];
		T eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
    		// decrease V
    		T estimate_error_Vx = 0;
    		{
				T eb_Vx_ = eb_Vx / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx_) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vx = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
            // decrease V
    		T estimate_error_Vy = 0;
    		{
				T eb_Vy_ = eb_Vy / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy_) + compute_bound_x_square(Vz[i], eb_Vz);
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vy = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
            // decrease V
    		T estimate_error_Vz = 0;
    		{
				T eb_Vz_ = eb_Vz / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz_);
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vz = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
    		// decrease P D
    		T estimate_error_P = 0;
    		{
				T eb_P_ = eb_P / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P_, eb_D);
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_P = compute_bound_division(V_TOT, C, e_V_TOT, e_C);    			
    		}
            T estimate_error_D = 0;
    		{
				T eb_D_ = eb_D / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D_);
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_D = compute_bound_division(V_TOT, C, e_V_TOT, e_C);    			
    		}
    		// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << " " << estimate_error_P << " " << estimate_error_D << std::endl;
			const T relative_epsilon = 1e-3;
            T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz, estimate_error_P, estimate_error_D});
            T epsilon = std::max(relative_epsilon * min_error, static_cast<T>(1e-12));
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            bool close_P  = fabs(estimate_error_P - min_error) < epsilon;
            bool close_D  = fabs(estimate_error_D - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
            if (close_P)  eb_P /= 1.5;
            if (close_D)  eb_D /= 1.5;
            if ((ebs[0] / eb_Vx > 10) || (ebs[1] / eb_Vy > 10) || (ebs[2] / eb_Vz > 10) || (ebs[3] / eb_P > 10) || (ebs[4] / eb_D > 10)) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}
    return true;
}

template<class T>
bool halfing_error_Mach_uniform(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T eb_P = ebs[3];
	T eb_D = ebs[4];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	T max_value = 0;
	int max_index = 0;
    T max_Vx = 0;
    T max_Vy = 0;
    T max_Vz = 0;
    T max_P = 0;
    T max_D = 0;
    T max_e_VTOT_2 = 0;
    T max_VTOT_2 = 0;
    T max_e_T = 0;
    T max_T = 0;
    T max_C = 0;
    T max_e_C = 0;
    T max_e_Mach = 0;
    T max_Mach = 0;
    int max_weight_Vx = 0;
    int max_weight_Vy = 0;
    int max_weight_Vz = 0;
    int max_weight_P = 0;
    int max_weight_D = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);
		// error of C
		T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		T C = c_2 * sqrt(Temp);
		// print_error("C", C, C_ori[i], e_C);
		T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		T Mach = V_TOT / C;
		// print_error("Mach", Mach, Mach_ori[i], e_Mach);
		error_est_Mach[i] = e_Mach;
		error_Mach[i] = Mach - Mach_ori[i];
		if(max_value < error_est_Mach[i]){
			max_value = error_est_Mach[i];
			max_index = i;
            max_Vx = Vx[i];
            max_Vy = Vy[i];
            max_Vz = Vz[i];
            max_P = P[i];
            max_D = D[i];
            max_e_VTOT_2 = e_V_TOT_2;
            max_VTOT_2 = V_TOT_2;
            max_e_T = e_T;
            max_T = Temp;
            max_e_C = e_C;
            max_C = C;
            max_e_Mach = e_Mach;
            max_Mach = Mach;
            max_weight_Vx = weights[0][i];
            max_weight_Vy = weights[1][i];
            max_weight_Vz = weights[2][i];
            max_weight_P = weights[3][i];
            max_weight_D = weights[4][i];
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
    // std::cout << "e_Mach = " << max_e_Mach << ", Mach = " << max_Mach << std::endl;
    // std::cout << "e_C = " << max_e_C << ", C = " << max_C << std::endl;
    // std::cout << "e_VTOT_2 = " << max_e_VTOT_2 << ", VTOT_2 = " << max_VTOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
    // std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;
    // std::cout << "e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    // std::cout << "max_weight_P = " << max_weight_P << ", max_weight_D = " << max_weight_D << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		T eb_P = ebs[3];
		T eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			T e_V_TOT_2 = 0;
			if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
			T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			T e_V_TOT = 0;
			if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			T V_TOT = sqrt(V_TOT_2);
			T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
			T Temp = P[i] / (D[i] * R);
			T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
			T C = c_2 * sqrt(Temp);
			estimate_error = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
            if (ebs[0] / eb_Vx > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}
    return true;
}

template<class T>
bool halfing_error_Mach_coordinate(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T eb_P = ebs[3];
	T eb_D = ebs[4];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	T max_value = 0;
	int max_index = 0;
    T max_Vx = 0;
    T max_Vy = 0;
    T max_Vz = 0;
    T max_P = 0;
    T max_D = 0;
    T max_e_VTOT_2 = 0;
    T max_VTOT_2 = 0;
    T max_e_T = 0;
    T max_T = 0;
    T max_C = 0;
    T max_e_C = 0;
    T max_e_Mach = 0;
    T max_Mach = 0;
    int max_weight_Vx = 0;
    int max_weight_Vy = 0;
    int max_weight_Vz = 0;
    int max_weight_P = 0;
    int max_weight_D = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);
		// error of C
		T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		T C = c_2 * sqrt(Temp);
		// print_error("C", C, C_ori[i], e_C);
		T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		T Mach = V_TOT / C;
		// print_error("Mach", Mach, Mach_ori[i], e_Mach);
		error_est_Mach[i] = e_Mach;
		error_Mach[i] = Mach - Mach_ori[i];
		if(max_value < error_est_Mach[i]){
			max_value = error_est_Mach[i];
			max_index = i;
            max_Vx = Vx[i];
            max_Vy = Vy[i];
            max_Vz = Vz[i];
            max_P = P[i];
            max_D = D[i];
            max_e_VTOT_2 = e_V_TOT_2;
            max_VTOT_2 = V_TOT_2;
            max_e_T = e_T;
            max_T = Temp;
            max_e_C = e_C;
            max_C = C;
            max_e_Mach = e_Mach;
            max_Mach = Mach;
            max_weight_Vx = weights[0][i];
            max_weight_Vy = weights[1][i];
            max_weight_Vz = weights[2][i];
            max_weight_P = weights[3][i];
            max_weight_D = weights[4][i];
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
    // std::cout << "e_Mach = " << max_e_Mach << ", Mach = " << max_Mach << std::endl;
    // std::cout << "e_C = " << max_e_C << ", C = " << max_C << std::endl;
    // std::cout << "e_VTOT_2 = " << max_e_VTOT_2 << ", VTOT_2 = " << max_VTOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
    // std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;
    // std::cout << "e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    // std::cout << "max_weight_P = " << max_weight_P << ", max_weight_D = " << max_weight_D << std::endl;

    // ************************************************ Vx, Vy, Vz, P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		T eb_P = ebs[3];
		T eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
    		// decrease V
    		T estimate_error_Vx = 0;
    		{
				T eb_Vx_ = eb_Vx / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx_ / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vx = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
            // decrease V
    		T estimate_error_Vy = 0;
    		{
				T eb_Vy_ = eb_Vy / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy_ / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vy = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
            // decrease V
    		T estimate_error_Vz = 0;
    		{
				T eb_Vz_ = eb_Vz / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz_ / static_cast<T>(std::pow(2.0, weights[2][i])));
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_Vz = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
    		}
    		// decrease P D
    		T estimate_error_P = 0;
    		{
				T eb_P_ = eb_P / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P_ / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_P = compute_bound_division(V_TOT, C, e_V_TOT, e_C);    			
    		}
            T estimate_error_D = 0;
    		{
				T eb_D_ = eb_D / 1.5;
				T e_V_TOT_2 = 0;
				if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
				T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
				T e_V_TOT = 0;
				if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				T V_TOT = sqrt(V_TOT_2);
				T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D_ / static_cast<T>(std::pow(2.0, weights[4][i])));
				T Temp = P[i] / (D[i] * R);
				T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
				T C = c_2 * sqrt(Temp);
				estimate_error_D = compute_bound_division(V_TOT, C, e_V_TOT, e_C);    			
    		}
    		// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << " " << estimate_error_P << " " << estimate_error_D << std::endl;
			const T relative_epsilon = 1e-3;
            T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz, estimate_error_P, estimate_error_D});
            T epsilon = std::max(relative_epsilon * min_error, static_cast<T>(1e-12));
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            bool close_P  = fabs(estimate_error_P - min_error) < epsilon;
            bool close_D  = fabs(estimate_error_D - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
            if (close_P)  eb_P /= 1.5;
            if (close_D)  eb_D /= 1.5;
            if ((ebs[0] / eb_Vx > 10) || (ebs[1] / eb_Vy > 10) || (ebs[2] / eb_Vz > 10) || (ebs[3] / eb_P > 10) || (ebs[4] / eb_D > 10)) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}
    return true;
}

template<class T>
T compute_global_value_range(const std::vector<T> data_vec){
	T global_max = 0, global_min = 0;
	T local_max = -std::numeric_limits<T>::max();
	T local_min = std::numeric_limits<T>::max();
	for(int i=0; i<data_vec.size(); i++){
		if(data_vec[i] > local_max) local_max = data_vec[i];
		if(data_vec[i] < local_min)	local_min = data_vec[i];
	}
	if(std::is_same<T, double>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	else if(std::is_same<T, float>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
	}
	return global_max - global_min;
}

template<class T>
std::vector<size_t> retrieve_Mach_Dummy(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::DummyApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3) reconstructors.back().mask = mask;
            reconstructors.back().load_metadata();
        }
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::DummyApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3) reconstructors.back().mask = mask;
            reconstructors.back().load_metadata();
            reconstructors.back().load_weight();
            reconstructors.back().span_weight();
            weights[i] = reconstructors.back().get_int_weights();
        }
        weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[n_variable - 1].get_weight_file_size();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Mach.data()), error_est_Mach.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_Mach_SZ3(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::SZApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3) reconstructors.back().mask = mask;
            reconstructors.back().load_metadata();
        }
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            else tolerance_met = halfing_error_Mach_coordinate(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::SZApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3){
                reconstructors.back().mask = mask;
                if (i==0) reconstructors.back().fetch_weight = true;
                else reconstructors.back().copy_int_weights(weights[0]);
            }
            else{
                if (i==3) reconstructors.back().fetch_weight = true;
                else reconstructors.back().copy_int_weights(weights[3]);
            }
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
            int index = (i < 3) ? 0 : 3;
            ebs[i] *= static_cast<T>(std::pow(2.0, reconstructors[index].get_max_weight()));
        }    
        weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[3].get_weight_file_size();
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                int index = (i < 3) ? 0 : 3;
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[index].get_max_weight())), -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            else tolerance_met = halfing_error_Mach_coordinate(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Mach.data()), error_est_Mach.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_Mach_PMGARD(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 9;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = DirectInterleaver<T>();
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
        }
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                if(i < 3){
                    // reconstruct with mask
                    int index = 0;
                    for(int j=0; j<num_elements; j++){
                        if(mask[j]){
                            reconstructed_vars[i][j] = reconstructed_data[index ++];
                        }
                        else{
                            reconstructed_vars[i][j] = 0; 
                            index++;
                        }
                    }
                }
                else{
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
    }
    else{
        std::vector<MDR::WeightReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, DirectInterleaver<int>, WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 9;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = DirectInterleaver<T>();
            auto weight_interleaver = DirectInterleaver<int>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, weight_interleaver, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
            reconstructors.back().load_weight();
            reconstructors.back().span_weight();
            weights[i] = reconstructors.back().get_int_weights();
        }
        weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[n_variable - 1].get_weight_file_size();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                if(i < 3){
                    // reconstruct with mask
                    int index = 0;
                    for(int j=0; j<num_elements; j++){
                        if(mask[j]){
                            reconstructed_vars[i][j] = reconstructed_data[index ++];
                        }
                        else{
                            reconstructed_vars[i][j] = 0; 
                            index++;
                        }
                    }
                }
                else{
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Mach.data()), error_est_Mach.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_Mach_GE(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::GEApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::GEApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3) reconstructors.back().mask = mask;
            reconstructors.back().load_metadata();
        }
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                // if(i < 3){
                //     // reconstruct with mask
                //     int index = 0;
                //     for(int j=0; j<num_elements; j++){
                //         if(mask[j]){
                //             reconstructed_vars[i][j] = reconstructed_data[index ++];
                //         }
                //         else{
                //             reconstructed_vars[i][j] = 0; 
                //             index++;
                //         }
                //     }
                // }
                // else{
                //     memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                // }
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            else tolerance_met = halfing_error_Mach_coordinate(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    else{
        std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::GEApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor_GE<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i < 3){
                reconstructors.back().mask = mask;
                if (i==0) reconstructors.back().fetch_weight = true;
                else reconstructors.back().copy_int_weights(weights[0]);
            }
            else{
                if (i==3) reconstructors.back().fetch_weight = true;
                else reconstructors.back().copy_int_weights(weights[3]);
            }
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
            int index = (i < 3) ? 0 : 3;
            ebs[i] *= static_cast<T>(std::pow(2.0, reconstructors[index].get_max_weight()));
        }
        weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[3].get_weight_file_size();
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                int index = (i < 3) ? 0 : 3;
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[index].get_max_weight())), -1);
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                // if(i < 3){
                //     // reconstruct with mask
                //     int index = 0;
                //     for(int j=0; j<num_elements; j++){
                //         if(mask[j]){
                //             reconstructed_vars[i][j] = reconstructed_data[index ++];
                //         }
                //         else{
                //             reconstructed_vars[i][j] = 0; 
                //             index++;
                //         }
                //     }
                // }
                // else{
                //     memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                // }
            }
            Vx_dec = reconstructed_vars[0].data();
            Vy_dec = reconstructed_vars[1].data();
            Vz_dec = reconstructed_vars[2].data();
            P_dec = reconstructed_vars[3].data();
            D_dec = reconstructed_vars[4].data();
            // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
            // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
            // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Mach = std::vector<T>(num_elements);
            error_est_Mach = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_Mach_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            else tolerance_met = halfing_error_Mach_coordinate(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Mach.data()), error_est_Mach.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Mach);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Mach);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    return total_retrieved_size;
}

int main(int argc, char ** argv){
    if(argc != 7){
		std::cout << "usage: para_VTOT [Approximator] [weighted] [decrease_method] [eb] [input path] [output path]" << std::endl;
 	}
	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	std::ostringstream oss;
	oss << rank;

    using T = double;
	int argv_id = 1;
    int compressor = atoi(argv[argv_id++]);
    int weighted = atoi(argv[argv_id++]);
    int decrease_method = atoi(argv[argv_id++]);
    T target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
    std::string output_path = argv[argv_id++];

    data_prefix_path += oss.str();
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<T> ebs; 
    ebs.push_back(compute_global_value_range(Vx_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(Vy_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(Vz_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(D_ori)*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(Vx_ori.data(), Vx_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(Vy_ori.data(), Vy_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(Vz_ori.data(), Vz_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(P_ori.data(), P_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(D_ori.data(), D_ori.size())*target_rel_eb);
	int n_variable = ebs.size();

    std::vector<T> Mach(num_elements);
    compute_Mach(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), P_ori.data(), D_ori.data(), num_elements, Mach.data());

    double tau, global_tau;
	global_tau = compute_global_value_range(Mach) * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
	Mach_ori = Mach.data();

    std::string mask_file = rdata_file_prefix + "mask.bin";
    uint32_t mask_file_size = 0;
    auto mask = readmask(mask_file.c_str(), mask_file_size);
    T local_max_act_error = 0;
	T local_max_est_error = 0;
    size_t weight_file_size = 0;
    std::vector<size_t> total_retrieved_size;
    switch (compressor)
	{
	case Dummy:
		total_retrieved_size = retrieve_Mach_Dummy(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	case SZ3:
		total_retrieved_size = retrieve_Mach_SZ3(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	case PMGARD:
		total_retrieved_size = retrieve_Mach_PMGARD(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
    case GE:
        total_retrieved_size = retrieve_Mach_GE(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	default:
		break;
	}
	double global_elapsed_time = 0;
    MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;

    if(!rank) printf("requested_error = %.10f\n", global_tau);
	
    double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_est_error = %.10f\n", global_max_est_error);

	double global_max_act_error = 0;
	MPI_Reduce(&local_max_act_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_act_error = %.10f\n", global_max_act_error);

	unsigned long long local_total_size = static_cast<unsigned long long>(mask_file_size) + std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0ULL) + static_cast<unsigned long long>(weight_file_size);

    // for(int i=0; i<size; i++){
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	if(rank == i) {
	// 		std::cout << "rank = " << rank << " act_iter = " << iter << ", local_bitrate = " << 8 * local_total_size * 1.0 / (num_elements * n_variable) << std::endl;
	// 	}
	// }

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);

    MPI_Finalize();

    return 0;
}