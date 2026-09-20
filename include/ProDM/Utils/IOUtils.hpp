#ifndef PRODM_UTILS_IOUTILS_HPP
#define PRODM_UTILS_IOUTILS_HPP
// IOUtils.hpp - library-wide helpers for raw binary files (readfile, readfile_pointer, writefile).
// Printing and statistics helpers are in StatUtils.hpp.
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
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
}  // namespace ProDM

#endif
