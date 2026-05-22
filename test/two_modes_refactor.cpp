#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "MDR/Refactor/Refactor.hpp"
#include "MDR/Tuner/Tuner.hpp"

using namespace std;
bool negabinary = true;
bool greedy = false;
bool bfs = false;
bool greedy_bfs = false;
std::vector<int> coeff_interp_directions;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor_time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Tuner>
std::vector<uint32_t> tune_interp_order(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, int stride, int block_size, Tuner tuner){
    struct timespec start, end;
    int err = 0;
    cout << "Start Tuning Per Level or Per Layer interpolation" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    tuner.tune(data.data(), dims, target_level, num_bitplanes, stride, block_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Tuner_time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    std::vector<uint32_t> interp_order = tuner.get_best_interp_order();
    std::cout << "best interp order = ";
    for(int i=0; i<interp_order.size(); i++){
        std::cout << interp_order[i] << " ";
    }
    std::cout << std::endl;
    return interp_order;
}

template<class T, class Level_Decomposer, class Layer_Decomposer, class Encoder, class Compressor, class Level_ErrorEstimator, class Layer_ErrorEstimator, class Level_SizeInterpreter, class Layer_SizeInterpreter>
std::vector<uint32_t> tune_interp_order_launcher(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, uint8_t mode, Level_Decomposer level_decomposer, Layer_Decomposer layer_decomposer, Encoder encoder, Compressor compressor, Level_ErrorEstimator level_estimator, Layer_ErrorEstimator layer_estimator, Level_SizeInterpreter level_interpreter, Layer_SizeInterpreter layer_interpreter){
    auto tuner  = MDR::ProfilingSamplingTuner<T, Level_Decomposer, Layer_Decomposer, Encoder, Compressor, Level_SizeInterpreter, Layer_SizeInterpreter, Level_ErrorEstimator, Layer_ErrorEstimator>(level_decomposer, layer_decomposer, encoder, compressor, level_interpreter, layer_interpreter);
    tuner.negabinary = negabinary;
    tuner.mode = mode;
    size_t num_elements = 0;
    return tune_interp_order(data, dims, target_level, num_bitplanes, 7, 7, tuner);
}

template<class T>
std::vector<uint32_t> tune_interp_order_initiator(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, std::string prior_mode){
    int num_dims = dims.size();
    uint8_t mode = (!strcmp(prior_mode.c_str(), "-PSNR"))? 1 : 0;
    using T_stream = typename std::conditional<std::is_same<T, float>::value, uint32_t, uint64_t>::type;
    auto level_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
    auto layer_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    negabinary = true;
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    auto level_estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
    auto level_interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(level_estimator);
    auto layer_estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
    auto layer_interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(layer_estimator);
    return tune_interp_order_launcher<T>(data, dims, target_level, num_bitplanes, mode, level_decomposer, layer_decomposer, encoder, compressor, level_estimator, layer_estimator, level_interpreter, layer_interpreter);
}

// Per level
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void init_fuse_composed_refactor(std::vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::FuseComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    // refactor.print();
    refactor.negabinary = negabinary;
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

// Per Layer
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void init_composed_refactor_new(std::vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor_new<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    // refactor.print();
    refactor.negabinary = negabinary;
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

// // Per level + CP
// template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
// void init_ordered_cp_refactor(std::vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
//     auto refactor = MDR::OrderedCPRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
//     refactor.print();
//     refactor.negabinary = negabinary;
//     refactor.coeff_interp_directions = coeff_interp_directions;
//     // refactor.greedy = greedy;
//     // refactor.bfs = bfs;
//     refactor.greedy_bfs = true;
//     evaluate(data, dims, target_level, num_bitplanes, refactor);
// }

// Per level + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void init_ordered_cp_refactor(std::vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::PartialOrderedCPRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    // refactor.print();
    refactor.negabinary = negabinary;
    refactor.coeff_interp_directions = coeff_interp_directions;
    // refactor.greedy = greedy;
    // refactor.bfs = bfs;
    // refactor.greedy_bfs = true;
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

// Per Layer + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void init_cp_refactor_new(std::vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::CPRefactor_new<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    // refactor.print();
    refactor.negabinary = negabinary;
    refactor.coeff_interp_directions = coeff_interp_directions;
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

void overall_refactoring_initiator(std::string filename, std::string output_path, std::string data_type, const vector<uint32_t>& dims, int target_level, int num_bitplanes, std::string prior_mode, std::string encoder_option, std::string cp){
    std::string metadata_file = output_path + "/metadata.bin";
    std::vector<std::string> files;
    for(int i=0; i<=(target_level + 2)*dims.size(); i++){
        string filename = output_path + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    std::string data_file = output_path + "/data.bin";
    if(!strcmp(data_type.c_str(), "-f")){
        using T = float;
        using T_stream = uint32_t;
        if(num_bitplanes > 32){
            num_bitplanes = 32;
            std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
        }
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        std::vector<uint32_t> best_interp_order = tune_interp_order_initiator(data, dims, target_level, num_bitplanes, prior_mode);
        if (best_interp_order[0] + best_interp_order[1] + best_interp_order[2] == 0){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto collector = MDR::SquaredErrorCollector<T>();
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            } else { // no CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            }
        } else { // Per Layer
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
            decomposer.interp_order = best_interp_order;
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto collector = MDR::SquaredErrorCollector<T>();
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            } else { // no CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            }
        }
    } else if(!strcmp(data_type.c_str(), "-d")){
        using T = double;
        using T_stream = uint64_t;
        if(num_bitplanes > 64){
            num_bitplanes = 64;
            std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
        }
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        std::vector<uint32_t> best_interp_order = tune_interp_order_initiator(data, dims, target_level, num_bitplanes, prior_mode);
        if (best_interp_order[0] + best_interp_order[1] + best_interp_order[2] == 0){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto collector = MDR::SquaredErrorCollector<T>();
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_ordered_cp_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            } else { // no CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_fuse_composed_refactor(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            }
        } else { // Per Layer
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
            decomposer.interp_order = best_interp_order;
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto collector = MDR::SquaredErrorCollector<T>();
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_cp_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            } else { // no CP
                auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    negabinary = true;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    negabinary = false;
                    init_composed_refactor_new(data, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
                }
            }
        }
    } else {
        std::cerr << "Only two float type supported: -f or -d" << std::endl;
        exit(-1);
    }
}

void usage(char* cmd) {
    std::cout << "two_modes_refactor usage: " << cmd <<
                  " data_file -[dataType: f/d] target_level num_bitplanes num_dims dim1 dim2 ... dimn output_path -[encoder_option: Nega/XOR/PerBit] -[prior_mode: eb(default)/PSNR] -[CP_or_not: CP/no_CP] (coeff_interp_direction, default tune)"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 -d 4 60 3 256 384 384 /output/path -PerBit -eb -CP" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string data_type = string(argv[argv_id ++]);
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
    int min_dim = *std::min_element(dims.begin(), dims.end()) - 1;
    int max_level = log2(min_dim);
    if(target_level > max_level){
        target_level = max_level;
        std::cout << "Target level exceeds the min dimension, change to target_level = " << target_level << " for correctness" << std::endl;
    }
    string output_path = string(argv[argv_id++]);
    string encoder_option = string(argv[argv_id++]);
    string mode = string(argv[argv_id++]);
    string cp = string(argv[argv_id++]);
    // std::cout << "argv_id = " << argv_id << ", argc = " << argc << std::endl;
    
    if(argv_id < argc){
        coeff_interp_directions.resize(num_dims);
        for(int i=0; i<num_dims; i++){
            coeff_interp_directions[i] = atoi(argv[argv_id++]);
        }
    }

    overall_refactoring_initiator(filename, output_path, data_type, dims, target_level, num_bitplanes, mode, encoder_option, cp);
    
    return 0;
}