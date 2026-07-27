#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "MDR/Tuner/Tuner.hpp"

using namespace std;
bool negabinary = false;

template <class T, class Tuner>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, int stride, int block_size, Tuner tuner){
    struct timespec start, end;
    int err = 0;
    cout << "Start Tuning" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    tuner.tune(data.data(), dims, target_level, num_bitplanes, stride, block_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Tuner time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    std::vector<uint32_t> interp_order = tuner.get_best_interp_order();
    std::cout << "best interp order = ";
    for(int i=0; i<interp_order.size(); i++){
        std::cout << interp_order[i] << " ";
    }
    std::cout << std::endl;
}

template<class T, class Level_Decomposer, class Layer_Decomposer, class Encoder, class Compressor, class Level_ErrorEstimator, class Layer_ErrorEstimator, class Level_SizeInterpreter, class Layer_SizeInterpreter>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, int stride, int block_size, Level_Decomposer level_decomposer, Layer_Decomposer layer_decomposer, Encoder encoder, Compressor compressor, Level_ErrorEstimator level_estimator, Layer_ErrorEstimator layer_estimator, Level_SizeInterpreter level_interpreter, Layer_SizeInterpreter layer_interpreter){
    // auto tuner = MDR::CoeffNaiveSamplingTuner<T, Decomposer, Encoder, Compressor, SizeInterpreter, ErrorEstimator>(decomposer, encoder, compressor, interpreter);
    auto tuner  = MDR::ProfilingSamplingTuner<T, Level_Decomposer, Layer_Decomposer, Encoder, Compressor, Level_SizeInterpreter, Layer_SizeInterpreter, Level_ErrorEstimator, Layer_ErrorEstimator>(level_decomposer, layer_decomposer, encoder, compressor, level_interpreter, layer_interpreter);
    tuner.negabinary = negabinary;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, stride, block_size, tuner);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int target_level = atoi(argv[argv_id ++]);
    int num_bitplanes = atoi(argv[argv_id ++]);
    if(num_bitplanes % 2 == 1) {
        num_bitplanes += 1;
        std::cout << "Change to " << num_bitplanes + 1 << " bitplanes for simplicity of negabinary encoding" << std::endl;
    }
    int num_dims = atoi(argv[argv_id ++]);
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id ++]);
    }
    int stride = atoi(argv[argv_id ++]);
    int block_size = atoi(argv[argv_id ++]);
    
    // using T = float;
    // using T_stream = uint32_t;
    // if(num_bitplanes > 32){
    //     num_bitplanes = 32;
    //     std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
    // }
    using T = double;
    using T_stream = uint64_t;
    if(num_bitplanes > 64){
        num_bitplanes = 64;
        std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
    }

    auto level_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
    auto layer_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    // negabinary = true;
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    auto level_estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
    // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    auto level_interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(level_estimator);
    auto layer_estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
    auto layer_interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(layer_estimator);
    // auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);

    test<T>(filename, dims, target_level, num_bitplanes, stride, block_size, level_decomposer, layer_decomposer, encoder, compressor, level_estimator, layer_estimator, level_interpreter, layer_interpreter);
    return 0;
}