#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Refactor/MDR/Refactor.hpp"
#include "mpi.h"

using namespace std;
bool negabinary = true;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    // struct timespec start, end;
    // int err = 0;
    // cout << "Start refactoring" << endl;
    // err = clock_gettime(CLOCK_REALTIME, &start);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double refactor_time = -MPI_Wtime();
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    refactor_time += MPI_Wtime();
    double global_elapsed_time = 0;
    MPI_Reduce(&refactor_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_elapsed_time = " << global_elapsed_time << std::endl;
    // err = clock_gettime(CLOCK_REALTIME, &end);
    // cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int argv_id = 1;
    string filename_base = string(argv[argv_id ++]);
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
    string output_path_base = string(argv[argv_id++]);

    // int tmp_rank = rank;
    // int idx = tmp_rank / 64;
    // tmp_rank = tmp_rank % 64;
    // int idy = tmp_rank / 8;
    // tmp_rank = tmp_rank % 8;
    // int idz = tmp_rank;

    // string filename = filename_base + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + ".d64";
    // string output_path = output_path_base + "/" + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + "/";

    string filename = filename_base + to_string(rank) + ".d64";
    string output_path = output_path_base + "/" + to_string(rank) + "/";

    string metadata_file = output_path + "/metadata_mdr.bin";
    vector<string> files;
    for(int i=0; i<=target_level; i++){
        string filename = output_path + "/level_" + to_string(i) + "_mdr.bin";
        files.push_back(filename);
    }
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

    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    // negabinary = true;
    auto encoder = MDR::PerBitBPEncoder_old<T, T_stream>();
    negabinary = false;
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();
    auto collector = MDR::SquaredErrorCollector<T>();
    auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
    // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);

    test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
    MPI_Finalize();
    return 0;
}