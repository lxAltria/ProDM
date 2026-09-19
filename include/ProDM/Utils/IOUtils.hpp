#ifndef PRODM_UTILS_IOUTILS_HPP
#define PRODM_UTILS_IOUTILS_HPP
// IOUtils.hpp - library-wide helpers for raw binary files and reconstruction statistics.
// (Formerly in the MGARDx utilities; they have no relation to the multilevel decomposition.)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include "ProDM/Namespace.hpp"

namespace ProDM {

// Read a whole binary file of Type values; num receives the element count (empty vector if the file is missing).
template <typename Type>
std::vector<Type> readfile(const char* file, size_t& num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        num = 0;
        return std::vector<Type>();
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    std::vector<Type> data(num_elements);
    fin.read(reinterpret_cast<char*>(data.data()), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}
// Same, into a malloc'ed buffer the caller frees (NULL if the file is missing).
template <typename Type>
Type* readfile_pointer(const char* file, size_t& num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        num = 0;
        return NULL;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    Type* data = (Type*)malloc(num_elements * sizeof(Type));
    fin.read(reinterpret_cast<char*>(data), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}
template <typename Type>
void writefile(const char* file, Type* data, size_t num_elements) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char*>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}
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

}  // namespace ProDM

#endif
