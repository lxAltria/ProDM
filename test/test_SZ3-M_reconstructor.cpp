#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "SZ3/api/sz.hpp"

using namespace std;

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
}

inline int find_index(double target_rel_eb, double& rel_eb){
    int i = 0;
    while(target_rel_eb < rel_eb){
        i ++;
        rel_eb /= 10;
    }
    return i;
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

template<class T>
size_t PSZ3_delta_reconstructor(string refactor_dict, const vector<T>& data, double tolerance, std::vector<T>& reconstructed_data){
    T value_range = compute_value_range(data);
    double file_eb = 0.1;
    int current_ind = -1;
    auto file_ind = find_index(tolerance / value_range, file_eb);
    size_t retrieved_size = 0;
    T * tmp_reconstructed_data = (T *) malloc(data.size() * sizeof(T));
    string filename = refactor_dict + "/SZ3_M_eb_" + std::to_string(file_ind) + ".bin";
    size_t n = 0;
    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
    retrieved_size += n;
    SZ3_decompress(cmpData.data(), n, tmp_reconstructed_data);
    memcpy(reconstructed_data.data(), tmp_reconstructed_data, sizeof(T) * data.size());
    return retrieved_size;
}


template <class T>
void evaluate(string refactor_dict, const vector<T>& data, const std::vector<double>& tolerance){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    T value_range = compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        std::vector<T> reconstructed_data(data.size(), 0);
        size_t retrieved_size = 0;
        cout << "Requested tolerance #" << i << " = " << tolerance[i] * value_range << endl;
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        retrieved_size = PSZ3_delta_reconstructor(refactor_dict, data, tolerance[i] * value_range, reconstructed_data);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        cout << "Retrieved data size = " << retrieved_size << endl;
        MGARD::print_statistics(data.data(), reconstructed_data.data(), data.size());
        cout << "Bitrate = " << (retrieved_size * 8.0) / data.size() << std::endl;
        cout << endl;
    }
}

template <class T>
void evaluate(string refactor_dict, const vector<T>& data, const vector<uint32_t>& dims, int num_snapshot){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    PSZ3_delta_reconstructor(refactor_dict, data, dims, num_snapshot);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void test(string filename, string refactor_dict, const vector<double>& tolerance){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(refactor_dict, data, tolerance);
}

void usage(char* cmd) {
    std::cout << "SZ3-M usage: " << cmd <<
                  " data_file refactor_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 3 1e-1 1e-2 1e-3 -d" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string refactor_dict = string(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string dtype = string(argv[argv_id++]);
    using T_stream = uint32_t;
    if (strcmp(dtype.c_str(), "-f") == 0){
        using T = float;
        test<T>(filename, refactor_dict, tolerance);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        using T = double;
        test<T>(filename, refactor_dict, tolerance);
    }
    return 0;
}