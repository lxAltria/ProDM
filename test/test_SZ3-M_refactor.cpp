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
char * SZ3_compress(size_t num_elements, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(num_elements);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
char * SZ3_compress_3D(size_t num_elements, uint32_t n1, uint32_t n2, uint32_t n3, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(n1, n2, n3);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
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
void PSZ3_delta_refactor(string output_dict, const vector<T>& data, const vector<uint32_t>& dims, int num_snapshot){
    T value_range = compute_value_range(data);
    std::vector<double> rel_ebs;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb * value_range);
    }
    std::vector<T> data_buffer(data);
    std::vector<T> dec_data_buffer(data);
    if (dims.size() == 3){
        for(int i=0; i<num_snapshot; i++){
            string filename = output_dict + "/SZ3_M_eb_" + std::to_string(i) + ".bin";
            size_t compressed_size = 0;
            auto compressed_data = SZ3_compress_3D(data.size(), size_t(dims[0]), size_t(dims[1]), size_t(dims[2]), data_buffer.data(), rel_ebs[i], compressed_size);
            MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
            SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
            free(compressed_data);
        }
    } else if (dims.size() == 1){
        for(int i=0; i<num_snapshot; i++){
            string filename = output_dict + "/SZ3_M_eb_" + std::to_string(i) + ".bin";
            size_t compressed_size = 0;
            auto compressed_data = SZ3_compress(data.size(), data_buffer.data(), rel_ebs[i], compressed_size);
            MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
            SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
            free(compressed_data);
        }
    }
}

template <class T>
void evaluate(string output_dict, const vector<T>& data, const vector<uint32_t>& dims, int num_snapshot){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    PSZ3_delta_refactor(output_dict, data, dims, num_snapshot);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void test(string filename, string output_dict, const vector<uint32_t>& dims, int num_snapshot){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(output_dict, data, dims, num_snapshot);
}

void usage(char* cmd) {
    std::cout << "SZ3-M usage: " << cmd <<
                  " data_file output_dict num_snapshot num_dim dim0 .. dimn -[dataType: f/d]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 18 3 256 384 384 -d" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string output_dict = string(argv[argv_id++]);
    int num_snapshot = atoi(argv[argv_id ++]);
    
    int num_dims = atoi(argv[argv_id ++]);
    // if (num_dims != 3){
    //     perror("Only 3D data.\n");
    // }
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id ++]);
    }
    
    string dtype = string(argv[argv_id++]);
    using T_stream = uint32_t;
    if (strcmp(dtype.c_str(), "-f") == 0){
        using T = float;
        test<T>(filename, output_dict, dims, num_snapshot);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        using T = double;
        test<T>(filename, output_dict, dims, num_snapshot);
    }
    return 0;
}