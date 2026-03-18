#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "MDR/Refactor/Refactor.hpp"

using namespace std;
bool negabinary = true;
bool greedy = false;
bool bfs = false;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // refactor.print();
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class ErrorEstimator, class Writer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, ErrorEstimator estimator, Writer writer){
    auto refactor = MDR::OrderedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, ErrorEstimator, Writer>(decomposer, interleaver, encoder, compressor, collector, estimator, writer);
    refactor.negabinary = negabinary;
    refactor.greedy = greedy;
    refactor.bfs = bfs;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

int main(int argc, char ** argv){

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
    string output_path = string(argv[argv_id++]);
    int interpreter_option = 0;
    interpreter_option = atoi(argv[argv_id ++]);
    if(interpreter_option == 0){
        greedy = true;
    } else if (interpreter_option == 1){
        bfs = true;
    }
    int encoder_option = 0;
    encoder_option = atoi(argv[argv_id ++]);

    string metadata_file = output_path +  "/refactored_data/metadata.bin";
    string data_file = output_path + "/refactored_data/data.bin";
    
    if(!strcmp(data_type.c_str(), "-f")){
        using T = float;
        using T_stream = uint32_t;
        if(num_bitplanes > 32){
            num_bitplanes = 32;
            std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
        }
        // auto decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        // auto interleaver = MDR::SFCInterleaver<T>();
        // auto interleaver = MDR::BlockedInterleaver<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
        // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        // negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
        // negabinary = false;
        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();
        auto collector = MDR::SquaredErrorCollector<T>();
        auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
        auto writer = MDR::OrderedFileWriter(metadata_file, data_file);

        if(encoder_option == 0){
            auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
            negabinary = true;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        } else if (encoder_option == 1){
            auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
            negabinary = true;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        } else {
            auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            negabinary = false;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        }
    } else if (!strcmp(data_type.c_str(), "-d")){
        using T = double;
        using T_stream = uint64_t;
        if(num_bitplanes > 64){
            num_bitplanes = 64;
            std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
        }
        // auto decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        // auto interleaver = MDR::SFCInterleaver<T>();
        // auto interleaver = MDR::BlockedInterleaver<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
        // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        // negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
        // negabinary = false;
        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();
        auto collector = MDR::SquaredErrorCollector<T>();
        auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
        auto writer = MDR::OrderedFileWriter(metadata_file, data_file);

        if(encoder_option == 0){
            auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
            negabinary = true;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        } else if (encoder_option == 1){
            auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
            negabinary = true;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        } else {
            auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            negabinary = false;
            test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, estimator, writer);
        }
    } else {
        std::cerr << "Only two float type supported: -f or -d" << std::endl;
        exit(-1);
    }
    return 0;
}