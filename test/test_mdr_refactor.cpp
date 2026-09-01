#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <sys/stat.h>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Refactor/MDR/Refactor.hpp"

using namespace std;
bool negabinary = true;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

template <class T, class T_stream>
void launch_refactor(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, int encoder_option, string metadata_file, const vector<string>& files){
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();
    auto collector = MDR::SquaredErrorCollector<T>();
    auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
    // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);

    if(encoder_option == 0){
        auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        negabinary = true;
        test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
    } else {
        auto encoder = MDR::PerBitBPEncoder_old<T, T_stream>();
        negabinary = false;
        test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
    }
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_file output_dict target_level num_bitplanes num_dim dim0 .. dimn [Encoder: NegaBinary-0, PerBit-1] -[dataType: f/d]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 4 60 3 256 384 384 0 -d" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string output_path = string(argv[argv_id++]);
    mkdir(output_path.c_str(), 0777);
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

    string metadata_file = output_path + "/metadata.bin";
    vector<string> files;
    for(int i=0; i<=target_level; i++){
        string filename = output_path + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

    int encoder_option = atoi(argv[argv_id++]);
    std::string dtype = string(argv[argv_id++]);
    if (strcmp(dtype.c_str(), "-f") == 0){
        if(num_bitplanes > 32){
            num_bitplanes = 32;
            std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
        }
        launch_refactor<float, uint32_t>(filename, dims, target_level, num_bitplanes, encoder_option, metadata_file, files);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        if(num_bitplanes > 64){
            num_bitplanes = 64;
            std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
        }
        launch_refactor<double, uint64_t>(filename, dims, target_level, num_bitplanes, encoder_option, metadata_file, files);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }
    return 0;
}
