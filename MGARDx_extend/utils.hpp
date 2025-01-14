#ifndef _MGARD_UTILS_HPP
#define _MGARD_UTILS_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>

namespace MGARD{

using namespace std;

template<typename Type>
std::vector<Type> readfile(const char *file, size_t &num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        return std::vector<Type>();
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    auto data = std::vector<Type>(num_elements);
    fin.read(reinterpret_cast<char *>(&data[0]), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}
template<typename Type>
Type * readfile_pointer(const char *file, size_t &num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        return NULL;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    Type * data = (Type *) malloc(num_elements * sizeof(Type));
    fin.read(reinterpret_cast<char *>(data), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}
template<typename Type>
void writefile(const char *file, Type *data, size_t num_elements) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}
template <class T>
void print(T * data, size_t n1, size_t n2, string s){
    cout << "Print data: " << s << endl;
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            cout << data[i * n2 + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
template <class T>
void print(const vector<vector<T>>& data){
    for(int i=0; i<data.size(); i++){
        for(int j=0; j<data[i].size(); j++){
            cout << j << ":" << data[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
template <class T>
void print_statistics(const T * data_ori, const T * data_dec, size_t data_size){
    double max_val = data_ori[0];
    double min_val = data_ori[0];
    double max_abs = fabs(data_ori[0]);
    for(int i=0; i<data_size; i++){
        if(data_ori[i] > max_val) max_val = data_ori[i];
        if(data_ori[i] < min_val) min_val = data_ori[i];
        if(fabs(data_ori[i]) > max_abs) max_abs = fabs(data_ori[i]);
    }
    double max_err = 0;
    int pos = 0;
    double mse = 0;
    for(int i=0; i<data_size; i++){
        double err = data_ori[i] - data_dec[i];
        mse += err * err;
        if(fabs(err) > max_err){
            pos = i;
            max_err = fabs(err);
        }
    }
    mse /= data_size;
    double psnr = 20 * log10((max_val - min_val) / sqrt(mse));
    cout << "Max value = " << max_val << ", min value = " << min_val << endl;
    cout << "Max error = " << max_err << ", pos = " << pos << endl;
    cout << "MSE = " << mse << ", PSNR = " << psnr << endl;
}
template <class T>
void print_statistics(const T * data_ori, const T * data_dec, size_t data_size, size_t compressed_size){
    print_statistics(data_ori, data_dec, data_size);
    cout << "Compression ratio = " << data_size * sizeof(T) * 1.0 / compressed_size << endl;
}
// compute dimensions for each level
/*
@params dims: dimensions
@params target_level: number of levels to perform
*/
vector<vector<size_t>> init_levels(const vector<size_t>& dims, size_t target_level){
    vector<vector<size_t>> level_dims;
    // compute n_nodal in each level
    for(int i=0; i<=target_level; i++){
        level_dims.push_back(vector<size_t>(dims.size()));
    }
    for(int i=0; i<dims.size(); i++){
        int n = dims[i];
        for(int j=0; j<=target_level; j++){
            level_dims[target_level - j][i] = n;
            n = (n >> 1) + 1;
        }
    }
    return level_dims;
}
template <class T>
inline T interp_cubic(T a, T b, T c, T d) {
    return (-a + 9 * b + 9 * c - d) / 16;
}
template <class T>
inline T interp_quad_1(T a, T b, T c) {
    return (3 * a + 6 * b - c) / 8;
}
template <class T>
inline T interp_quad_2(T a, T b, T c) {
    return (-a + 6 * b + 3 * c) / 8;
}
}
#endif