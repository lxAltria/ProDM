#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "MDR/Refactor/Refactor.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "MDR/Interleaver/NewInterleaver.hpp"
#include "MDR/RefactorUtils.hpp"
#include "MDR/BitplaneEncoder/NegaBinaryBPEncoder.hpp"
#include "MDR/LosslessCompressor/LevelCompressor.hpp"

using namespace std;

template <class T, class T_stream>
void test_interleave(vector<T>& data, vector<vector<T>> & buffer, std::vector<std::vector<uint8_t *>> & streams_, std::vector<T> & level_error_bounds, const vector<uint32_t>& dims, int target_level){
    MDR::MGARDHierarchicalDecomposer_new<T> decomposer;
    decomposer.decompose(data.data(), dims, target_level);
    std::vector<std::vector<uint32_t>> level_dims = MDR::compute_level_dims_new(dims, target_level);
    std::vector<uint32_t> level_elements = MDR::compute_level_elements(level_dims, target_level);
    auto interleaver = MDR::DirectInterleaver_new<T>();
    std::vector<uint32_t> dims_dummy(dims.size(), 0);
    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    for(int i=0; i<=target_level; i++){
        const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        interleaver.interleave(data.data(), dims, level_dims[i], prev_dims, buffer[i].data(), i, target_level);
        T level_max_error = MDR::compute_max_abs_value(buffer[i].data(), level_elements[i]);
        level_error_bounds.push_back(level_max_error * 4);
        int level_exp = 0;
        frexp(level_max_error, &level_exp);
        std::vector<uint32_t> stream_sizes;
        auto streams = encoder.encode(buffer[i].data(), level_elements[i], level_exp, 16, stream_sizes);
        streams_.push_back(streams);
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Interleave time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class T_stream>
void test_reposistion(vector<T>& data, vector<vector<T>> & buffer, std::vector<std::vector<uint8_t *>> & streams_, std::vector<T> & level_error_bounds, const vector<uint32_t>& dims, int target_level){
    std::vector<std::vector<uint32_t>> level_dims = MDR::compute_level_dims_new(dims, target_level);
    std::vector<uint32_t> level_elements = MDR::compute_level_elements(level_dims, target_level);
    auto interleaver = MDR::DirectInterleaver_new<T>();
    std::vector<uint32_t> dims_dummy(dims.size(), 0);
    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    std::vector<uint8_t> BP = {17, 11, 10};
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    for(int i=0; i<=target_level; i++){
        int level_exp = 0;
        frexp(level_error_bounds[i] / 4, &level_exp);
        std::vector<const uint8_t *> const_streams(streams_[i].begin(), streams_[i].end());
        auto level_decoded_data = encoder.progressive_decode(const_streams, level_elements[i], level_exp, 0, 16, i);
        const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        // interleaver.reposition(buffer[i].data(), dims, level_dims[i], prev_dims, data.data(), i, target_level);
        interleaver.reposition(level_decoded_data, dims, level_dims[i], prev_dims, data.data(), i, target_level);
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Recomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    MDR::MGARDHierarchicalDecomposer_new<T> recomposer;
    recomposer.recompose(data.data(), dims, target_level);
}

template <class T>
int print_statistics(const T * data_ori, const T * data_dec, size_t data_size){
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
    return pos;
}

template <class T, class T_stream>
void test(string filename, const vector<uint32_t>& dims, int target_level){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    auto data_ori(data);
    vector<vector<T>> buffer;
    auto level_dims = MDR::compute_level_dims_new(dims, target_level);
    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
    for(int i=0; i<level_elements.size(); i++){
        buffer.push_back(vector<T>(level_elements[i], 0));
    }
    std::vector<std::vector<uint8_t *>> streams;
    std::vector<T> level_error_bounds;
    test_interleave<T, T_stream>(data, buffer, streams, level_error_bounds, dims, target_level);
    std::vector<T> repositioned_data(num_elements, 0);
    test_reposistion<T, T_stream>(repositioned_data, buffer, streams, level_error_bounds, dims, target_level);
    auto pos = print_statistics(data_ori.data(), repositioned_data.data(), num_elements);
    std::cout << data_ori[pos] << " " << repositioned_data[pos] << std::endl;
    std::cout << pos / (dims[1]*dims[2]) << " ";
    pos = pos % (dims[1]*dims[2]);
    std::cout << pos / dims[2] << " ";
    pos = pos % dims[2];
    std::cout << pos << std::endl;
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    int target_level = atoi(argv[3]);
    const int num_dims = atoi(argv[4]);
    vector<uint32_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[5 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    switch(type){
        case 0:
            {
                test<float, uint32_t>(filename, dims, target_level);
                break;
            }
        case 1:
            {
                test<double, uint64_t>(filename, dims, target_level);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}