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

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, const int block_size, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, target_level, num_bitplanes, block_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Decomposer, class InterleaverT, class InterleaverInt, class Encoder, class Compressor, class ErrorCollector, class Writer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, const int block_size, Decomposer decomposer, InterleaverT interleaver, InterleaverInt weight_interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::WeightRefactor<T, Decomposer, InterleaverT, InterleaverInt, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, weight_interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, block_size, refactor);
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
    const int block_size = atoi(argv[argv_id ++]);

    string metadata_file = "refactored_weight_cubic_data/metadata.bin";
    vector<string> files;
    for(int i=0; i<=target_level; i++){
        string filename = "refactored_weight_cubic_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    using T = float;
    using T_stream = uint32_t;
    if(num_bitplanes > 32){
        num_bitplanes = 32;
        std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
    }
    // using T = double;
    // using T_stream = uint64_t;
    // if(num_bitplanes > 64){
    //     num_bitplanes = 64;
    //     std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
    // }

    auto decomposer = MDR::MGARDCubicDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto weight_interleaver = MDR::DirectInterleaver<int>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    auto encoder = MDR::WeightedPerBitBPEncoder<T, T_stream>();
    negabinary = true;
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    // negabinary = false;
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();
    auto collector = MDR::SquaredErrorCollector<T>();
    auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
    // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);

    test<T>(filename, dims, target_level, num_bitplanes, block_size, decomposer, interleaver, weight_interleaver, encoder, compressor, collector, writer);
    return 0;
}