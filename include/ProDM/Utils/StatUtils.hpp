#ifndef PRODM_UTILS_STATUTILS_HPP
#define PRODM_UTILS_STATUTILS_HPP
// StatUtils.hpp - library-wide printing and statistics helpers: value ranges, reconstruction statistics
// (max error, MSE, PSNR, compression ratio), small vector printers, and the interpolation-coefficient normalizer.
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include "ProDM/Namespace.hpp"

namespace ProDM {

template <class T>
void print(T* data, size_t n1, size_t n2, std::string s) {
    std::cout << "Print data: " << s << std::endl;
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) std::cout << data[i * n2 + j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
template <class T>
void print(const std::vector<std::vector<T>>& data) {
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data[i].size(); j++) std::cout << j << ":" << data[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
// Value range, maximal error and its position, MSE, PSNR and NRMSE of a reconstruction against the original.
template <class T>
void print_statistics(const T* data_ori, const T* data_dec, size_t data_size) {
    double max_val = data_ori[0], min_val = data_ori[0], max_abs = std::fabs(data_ori[0]);
    for (size_t i = 0; i < data_size; i++) {
        if (data_ori[i] > max_val) max_val = data_ori[i];
        if (data_ori[i] < min_val) min_val = data_ori[i];
        if (std::fabs(data_ori[i]) > max_abs) max_abs = std::fabs(data_ori[i]);
    }
    double max_err = 0, mse = 0; size_t pos = 0;
    for (size_t i = 0; i < data_size; i++) {
        double err = data_ori[i] - data_dec[i];
        mse += err * err;
        if (std::fabs(err) > max_err) { pos = i; max_err = std::fabs(err); }
    }
    mse /= data_size;
    double psnr = 20 * std::log10((max_val - min_val) / std::sqrt(mse));
    std::cout << "Max value = " << max_val << ", min value = " << min_val << std::endl;
    std::cout << "Max error = " << max_err << ", pos = " << pos << std::endl;
    std::cout << "MSE = " << mse << ", PSNR = " << psnr << ", NRMSE = " << std::sqrt(mse) / (max_val - min_val) << std::endl;
}
template <class T>
void print_statistics(const T* data_ori, const T* data_dec, size_t data_size, size_t compressed_size) {
    print_statistics(data_ori, data_dec, data_size);
    std::cout << "Compression ratio = " << data_size * sizeof(T) * 1.0 / compressed_size << std::endl;
}

template <class T>
void print_error(std::string varname, T dec, T ori, T est){
	std::cout << varname << ": dec = " << dec << ", ori = " << ori << ", error = " << dec - ori << ", est = " << est << std::endl; 
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
T compute_value_range(const T* vec, uint32_t size){
	T min = vec[0];
	T max = vec[0];
	for(int i=0; i<size; i++){
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

}  // namespace ProDM

#endif
