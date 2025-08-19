#ifndef _MDR_QOI_UTILS_HPP
#define _MDR_QOI_UTILS_HPP

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <cstdint>

namespace MDR{

const std::vector<std::string> names{"V_TOT", "T", "C", "Mach", "PT", "mu"};

// f(x) = x^2
template <class T>
inline double compute_bound_x_square(T x, T eb){
	return 2*fabs(x)*eb + eb*eb;
}

// f(x) = x^2
template <class T>
inline double compute_inverse_bound_x_square(T x, T eb, T tau){
	return sqrt(x*x + tau) - x*x;
}

// f(x) = sqrt(x)
// template <class T>
// inline double compute_bound_square_root_x(T x, T eb){
// 	if(x == 0) {
// 		return sqrt(eb);
// 	}
// 	if(x > eb){
// 		return eb / (sqrt(x - eb) + sqrt(x));
// 	}
// 	else{
// 		// return eb / sqrt(x);
// 		T tau_1 = sqrt(x + eb) - sqrt(x);
// 		T tau_2 = sqrt(x);
// 		return (tau_1 > tau_2) ? tau_1 : tau_2;
// 	}
// }
template <class T>
inline double compute_bound_square_root_x(T x, T eb){
	if(x == 0) {
		return sqrt(eb);
	}
	if(x > eb){
		return eb / (sqrt(x - eb) + sqrt(x));
	}
	else{
		return eb / sqrt(x);
	}
}

// f(x) = sqrt(x)
template <class T>
inline double compute_inverse_bound_square_root_x(T x, T eb, T tau){
	if(x == 0){
		return tau * tau;
	}
	if(x > eb){
		return tau * (sqrt(x - eb) + sqrt(x));
	}
	else{
		return tau * sqrt(x);
	}
}

// f(x) = 1/(x+a)
template <class T>
inline double compute_bound_radical(T x, T a, T eb){
	if(fabs(x+a) > eb){
		if(x+a > 0){
			return eb / ((x+a-eb)*fabs(x+a));
		}
		else{
			return eb / ((x+a+eb)*fabs(x+a));
		}
	}
	else{
		std::cout << "Warning: cannot control error in 1/(x+a)\n";
		return 0;		
	}
}

// f(x) = 1/(x+a)
template <class T>
inline double compute_inverse_bound_radical(T x, T a, T eb, T tau){
	if(fabs(x+a) > eb){
		if(x+a > 0){
			return tau * ((x+a-eb)*fabs(x+a));
		}
		else{
			return tau * ((x+a+eb)*fabs(x+a));
		}
	}
	else{
		std::cout << "Warning: cannot control error in 1/(x+a)\n";
		return 0;		
	}
}

template <class T>
inline double compute_bound_multiplication(T x, T y, T eb_x, T eb_y){
	return fabs(x)*eb_y + fabs(y)*eb_x + eb_x*eb_y;
}

// f(x, y) = x/y
template <class T>
inline double compute_bound_division(T x, T y, T eb_x, T eb_y){
	if(eb_y < fabs(y)){
		double e = fabs(x)*eb_y + fabs(y)*eb_x;
		if(y > 0) return e / (y*(y - eb_y));
		else return e / (y*(y + eb_y));
	}
	else{
		// std::cout << "Warning: cannot control error in x/y\n";
		return 0;
	}
}

template <class T>
void print_error(std::string varname, T dec, T ori, T est){
	std::cout << varname << ": dec = " << dec << ", ori = " << ori << ", error = " << dec - ori << ", est = " << est << std::endl; 
}

template <class T>
double UniformTighteningVTOT(T Vx, T Vy, T Vz, T eb_Vx, T eb_Vy, T eb_Vz, T tau){
	// error of total velocity square
	double e_V_TOT_2 = compute_bound_x_square(Vx, eb_Vx) + compute_bound_x_square(Vy, eb_Vy) + compute_bound_x_square(Vz, eb_Vz);
	double V_TOT_2 = Vx*Vx + Vy*Vy + Vz*Vz;
	// error of total velocity
	double e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
	double V_TOT = sqrt(V_TOT_2);
	double new_e_V_TOT_2 = compute_inverse_bound_square_root_x(V_TOT, e_V_TOT, tau);
	// compute assignment on Vx, Vy, Vz
	double e1 = 2*fabs(Vx)*eb_Vx + 2*fabs(Vy)*eb_Vy + 2*fabs(Vz)*eb_Vz;
	double e2 = eb_Vx*eb_Vx + eb_Vy*eb_Vy  + eb_Vz*eb_Vz;
	// assuming the same scale on error bound deduction
	double scale = (-e1 + sqrt(e1*e1 + 4*e2*new_e_V_TOT_2))/(2*e2);
	return scale;
}

template <class T>
double UniformTighteningT(T P, T D, T eb_P, T eb_D, T c_1, T tau){
	double e_T = c_1 * compute_bound_division(P, D, eb_P, eb_D);
	return tau / e_T;	
}

template <class T>
void compute_QoIs(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n,
					T * V_TOT_, T * Temp_, T * C_, T * Mach_, T * PT_, T * mu_){
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	for(int i=0; i<n; i++){
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double Temp = P[i] / (D[i] * R);
		double C = sqrt(gamma * R * Temp);
		double Mach = V_TOT / C;
		double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
		double Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
		double PT = P[i] * Mach_tmp_mi;
		double TrS_TS = (T_r + S) / (Temp + S);
		double T_Tr_3 = pow(Temp/T_r, 3);
		double T_Tr_3_sqrt = sqrt(T_Tr_3);
		double mu = mu_r * T_Tr_3_sqrt * TrS_TS;
		V_TOT_[i] = V_TOT;
		Temp_[i] = Temp;
		C_[i] = C;
		Mach_[i] = Mach;
		PT_[i] = PT;
		mu_[i] = mu;
	}
}

template <class T>
void compute_VTOT(const T * Vx, const T * Vy, const T * Vz, size_t n, T * V_TOT_){
	for(int i=0; i<n; i++){
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		V_TOT_[i] = V_TOT;
	}
}

template <class T>
void compute_T(const T * P, const T * D, size_t n, T * Temp_){
	double R = 287.1;
	double gamma = 1.4;
	for(int i=0; i<n; i++){
		double Temp = P[i] / (D[i] * R);
		Temp_[i] = Temp;
	}
}

template <class T>
void compute_C(const T * P, const T * D, size_t n, T * C_){
	double R = 287.1;
	double gamma = 1.4;
	for(int i=0; i<n; i++){
		double Temp = P[i] / (D[i] * R);
		double C = sqrt(gamma * R * Temp);
		C_[i] = C;
	}
}

template <class T>
void compute_mu(const T * P, const T * D, size_t n, T * mu_){
	double R = 287.1;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	for(int i=0; i<n; i++){
		double Temp = P[i] / (D[i] * R);
		double TrS_TS = (T_r + S) / (Temp + S);
		double T_Tr_3 = pow(Temp/T_r, 3);
		double T_Tr_3_sqrt = sqrt(T_Tr_3);
		double mu = mu_r * T_Tr_3_sqrt * TrS_TS;
		mu_[i] = mu;
	}
}

template <class T>
void compute_Mach(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, T * Mach_){
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	for(int i=0; i<n; i++){
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double Temp = P[i] / (D[i] * R);
		double C = sqrt(gamma * R * Temp);
		double Mach = V_TOT / C;
		Mach_[i] = Mach;
	}
}

template <class T>
void compute_PT(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, T * PT_){
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	for(int i=0; i<n; i++){
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double Temp = P[i] / (D[i] * R);
		double C = sqrt(gamma * R * Temp);
		double Mach = V_TOT / C;
		double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
		double Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
		double PT = P[i] * Mach_tmp_mi;
		PT_[i] = PT;
	}
}

template <class T>
T compute_value_range(const std::vector<T>& vec){
	T min = vec[0];
	T max = vec[0];
	for(int i=0; i<vec.size(); i++){
		if(vec[i] < min) min = vec[i];
		if(vec[i] > max) max = vec[i];
	}
	return max - min;
}

template <class T>
void print_info(const std::string& name, const std::vector<T>& vec){
	T max = vec[0];
	T min = vec[0];
	for(int i=1; i<vec.size(); i++){
		if(max < vec[i]) max = vec[i];
		if(min > vec[i]) min = vec[i];
	}
	std::cout << name << ": min value = " << min << ", max value = " << max << ", value range = " << max - min << std::endl;
}

// template <class T>
// void print_max_abs(const std::string& name, const std::vector<T>& vec){
// 	T max = fabs(vec[0]);
// 	for(int i=1; i<vec.size(); i++){
// 		if(max < fabs(vec[i])) max = fabs(vec[i]);
// 	}
// 	std::cout << name << ": max absolute value = " << max << std::endl;
// }

template <class T>
T print_max_abs(const std::string& name, const std::vector<T>& vec){
	T max = fabs(vec[0]);
	for(int i=1; i<vec.size(); i++){
		if(max < fabs(vec[i])) max = fabs(vec[i]);
	}
	// std::cout << name << ": max absolute value = " << max << std::endl;
	return max;
}

template <class T>
inline void normalize_coefficient(T& coeff_0, T& coeff_1){
	double min_val = std::min(coeff_0, coeff_1);
	coeff_0 /= min_val;
	coeff_1 /= min_val;
	if(coeff_0 + coeff_1 == 2){
		coeff_0 = 1.5;
		coeff_1 = 1.5;
	}
}

}
#endif