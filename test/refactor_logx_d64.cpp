#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "PDR/Refactor/Refactor.hpp"
#define Dummy_Cmp 0
#define MGARD_Cmp 1
#define SZ2_Cmp 2
#define SZ3_Cmp 3
#define HPEZ_Cmp 4
#define GE_Cmp 5

using namespace std;
bool negabinary = true;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int num_bitplanes, Refactor refactor, int max_weight, int block_size, T approximator_eb){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, 0, num_bitplanes, approximator_eb, max_weight, block_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Approximator, class Encoder, class Compressor, class Writer>
void test(string filename, const vector<uint32_t>& dims, int num_bitplanes, Approximator approximator, Encoder encoder, Compressor compressor, Writer writer, int max_weight, int block_size, T approximator_eb){
    auto refactor = PDR::WeightedApproximationBasedRefactor<T, Approximator, Encoder, Compressor, Writer>(approximator, encoder, compressor, writer);
    refactor.negabinary = negabinary;
    refactor.store_weight = true;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::vector<T> weights(num_elements);
    for(int i=0; i<num_elements; i++){
        weights[i] = 1 / data[i];
    }
    refactor.QoI = weights;
    evaluate(data, dims, num_bitplanes, refactor, max_weight, block_size, approximator_eb);
}

template <class T, class Encoder, class Compressor, class Writer>
void launch_refactor(string filename, const vector<uint32_t>& dims, int num_bitplanes, int approximator_rank, Encoder encoder, Compressor compressor, Writer writer, int max_weight, int block_size, T approximator_eb){
    switch(approximator_rank){
        case Dummy_Cmp:{
            auto approximator = PDR::DummyApproximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        case MGARD_Cmp:{
            auto approximator = PDR::MGARDApproximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        case SZ2_Cmp:{
            auto approximator = PDR::SZ2Approximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        case SZ3_Cmp:{
            auto approximator = PDR::SZ3Approximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        case HPEZ_Cmp:{
            auto approximator = PDR::HPEZApproximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        case GE_Cmp:{
            auto approximator = PDR::GEApproximator<T>();
            test<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
            break;
        }
        default:
            perror("Undefined Approximator\n");
            break;
    }
}

void usage(char* cmd) {
    std::cout << "QPro usage: " << cmd <<
                  " data_file output_dict num_bitplanes num_dim dim0 .. dimn -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4] max_weight block_size approximator_eb"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 3 256 384 384 -d 4" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    std::string filename = string(argv[argv_id ++]);
    std::string output_dict = string(argv[argv_id++]);
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

    int target_level = 0; // #level = 1 for PDR
    std::string metadata_file = output_dict + "/metadata.bin";
    vector<string> files;
    for(int i=0; i<=target_level; i++){
        std::string filename = output_dict + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    
    using T_stream = uint32_t;

    std::string dtype = string(argv[argv_id++]);
    int approximator = atoi(argv[argv_id++]);
    int max_weight = atoi(argv[argv_id++]);
    int block_size = atoi(argv[argv_id++]);
    if (strcmp(dtype.c_str(), "-f") == 0){
        if(num_bitplanes > 32){
            num_bitplanes = 32;
            std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
        }
        using T = float;
        T approximator_eb = atoi(argv[argv_id++]);
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, T_stream>();
        negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        launch_refactor<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        if(num_bitplanes > 64){
            num_bitplanes = 64;
            std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
        }
        using T = double;
        T approximator_eb = atoi(argv[argv_id++]);
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, T_stream>();
        negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        launch_refactor<T>(filename, dims, num_bitplanes, approximator, encoder, compressor, writer, max_weight, block_size, approximator_eb);
    }

    return 0;
}