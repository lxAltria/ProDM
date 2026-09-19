#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>
#include <sys/stat.h>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Refactor/PDR/Refactor.hpp"
#include "ProDM/Utils/MaskUtils.hpp"
#define BP 0
#define WBP 1
#define QoI_Vtot 0
#define QoI_Vtot2 1

using namespace std;

const vector<string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};

template <class T>
void launch_refactor(const string& data_dir, const string& refactor_dir, int num_bitplanes, const vector<uint32_t>& dims, const string& suffix, int mode, int qoi, T approximator_eb, int max_weight, int block_size){
#ifdef PRODM_HAVE_HPEZ
    size_t num_elements = 1;
    for(const auto& dim:dims) num_elements *= dim;

    int n_variable = var_list.size();
    vector<vector<T>> vars_vec(n_variable);
    for(int i=0; i<n_variable; i++){
        string filename = data_dir + "/" + var_list[i] + suffix;
        size_t num = 0;
        vars_vec[i] = ProDM::readfile<T>(filename.c_str(), num);
        if(num != num_elements){
            cerr << "File " << filename << " has " << num << " elements, but the given dimensions require " << num_elements << endl;
            exit(-1);
        }
    }

    // mask of nonzero velocities
    vector<unsigned char> mask(num_elements, 0);
    for(size_t i=0; i<num_elements; i++){
        if(vars_vec[0][i]*vars_vec[0][i] + vars_vec[1][i]*vars_vec[1][i] + vars_vec[2][i]*vars_vec[2][i] != 0){
            mask[i] = 1;
        }
    }
    string mask_file = refactor_dir + "/mask.bin";
    ProDM::writemask(mask_file.c_str(), mask.data(), mask.size());

    uint8_t target_level = 0;

    if(mode == BP){
        for(int i=0; i<n_variable; i++){
            string rdir_prefix = refactor_dir + "/" + var_list[i] + "_refactored";
            mkdir(rdir_prefix.c_str(), 0777);
            string metadata_file = rdir_prefix + "/metadata.bin";
            vector<string> files = {rdir_prefix + "/level_0.bin"};

            auto approximator = PDR::HPEZApproximator<T>();
            auto encoder = ProDM::NegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = ProDM::AdaptiveLevelCompressor(64);
            auto writer = ProDM::ConcatLevelFileWriter(metadata_file, files);
            auto refactor = PDR::ApproximationBasedRefactor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(writer)>(approximator, encoder, compressor, writer);
            refactor.negabinary = true;
            refactor.mask = mask;
            refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, approximator_eb);
        }
    }
    else{
        // per-point weights derived from the requested QoI
        vector<T> weights(num_elements, 0);
        if(qoi == QoI_Vtot){
            // weight = 1 / Vtot
            for(size_t i=0; i<num_elements; i++){
                if(mask[i]){
                    T V = sqrt(vars_vec[0][i]*vars_vec[0][i] + vars_vec[1][i]*vars_vec[1][i] + vars_vec[2][i]*vars_vec[2][i]);
                    weights[i] = 1.0 / V;
                }
            }
        }
        else{
            // weight = |Vx| + |Vy| + |Vz|
            for(size_t i=0; i<num_elements; i++){
                T sum = 0;
                for(int j=0; j<n_variable; j++){
                    sum += fabs(vars_vec[j][i]);
                }
                weights[i] = sum;
            }
        }

        // the integer weights are derived once (first variable, which also stores
        // weight.bin in its directory) and shared by the other variables
        vector<int> int_weights;
        for(int i=0; i<n_variable; i++){
            string rdir_prefix = refactor_dir + "/" + var_list[i] + "_refactored";
            mkdir(rdir_prefix.c_str(), 0777);
            string metadata_file = rdir_prefix + "/metadata.bin";
            vector<string> files = {rdir_prefix + "/level_0.bin"};

            auto approximator = PDR::HPEZApproximator<T>();
            auto encoder = ProDM::WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = ProDM::AdaptiveLevelCompressor(64);
            auto writer = ProDM::ConcatLevelFileWriter(metadata_file, files);
            auto refactor = PDR::WeightedApproximationBasedRefactor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(writer)>(approximator, encoder, compressor, writer);
            refactor.negabinary = true;
            refactor.mask = mask;
            if(i == 0){
                refactor.QoI = weights;
                refactor.store_weight = true;
            }
            else{
                refactor.copy_int_weights(int_weights);
            }
            refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, approximator_eb, max_weight, block_size);
            if(i == 0) int_weights = refactor.get_int_weights();
        }
    }
#else
    cerr << "This tool requires the HPEZ approximator; rebuild with -DPRODM_WITH_HPEZ=ON" << endl;
    exit(-1);
#endif
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_dir refactor_dir num_bitplanes num_dim dim0 .. dimn -[dataType: f/d] [mode: BP-0, WBP-1] [QoI: Vtot-0, Vtot2-1] [approximator_eb max_weight block_size]"
                  << std::endl
                  << "reads VelocityX/Y/Z (.dat for -d, .dat.f32 for -f) from data_dir" << std::endl
                  << "example: " << cmd <<
                  " data refactor 60 3 100 500 500 -d 1 0 0.001 4 1" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    std::string data_dir = string(argv[argv_id++]);
    std::string refactor_dir = string(argv[argv_id++]);
    mkdir(refactor_dir.c_str(), 0777);
    int num_bitplanes = atoi(argv[argv_id++]);
    if(num_bitplanes % 2 == 1) {
        num_bitplanes += 1;
        std::cout << "Change to " << num_bitplanes << " bitplanes for simplicity of negabinary encoding" << std::endl;
    }
    int num_dims = atoi(argv[argv_id++]);
    if(num_dims <= 0 || argc < argv_id + num_dims + 3){
        std::cerr << "Insufficient or invalid arguments (num_dims parsed as " << num_dims << "); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id++]);
    }

    std::string dtype = string(argv[argv_id++]);
    int mode = atoi(argv[argv_id++]);
    int qoi = atoi(argv[argv_id++]);
    if((mode != BP && mode != WBP) || (qoi != QoI_Vtot && qoi != QoI_Vtot2)){
        std::cerr << "Unknown mode " << mode << " (expected BP-0 or WBP-1) or QoI " << qoi << " (expected Vtot-0 or Vtot2-1)" << std::endl;
        usage(argv[0]);
        return -1;
    }
    if(mode == WBP && num_dims > 3){
        std::cerr << "The weighted mode supports up to 3 dimensions" << std::endl;
        return -1;
    }
    double approximator_eb = 0.001;
    int max_weight = 4;
    int block_size = 1;
    if(argc > argv_id) approximator_eb = atof(argv[argv_id++]);
    if(argc > argv_id) max_weight = atoi(argv[argv_id++]);
    if(argc > argv_id) block_size = atoi(argv[argv_id++]);

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    if (strcmp(dtype.c_str(), "-f") == 0){
        if(num_bitplanes > 32){
            num_bitplanes = 32;
            std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
        }
        launch_refactor<float>(data_dir, refactor_dir, num_bitplanes, dims, ".dat.f32", mode, qoi, (float)approximator_eb, max_weight, block_size);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        if(num_bitplanes > 60){
            num_bitplanes = 60;
            std::cout << "Only less than 60 bitplanes are supported for double-precision floating point" << std::endl;
        }
        launch_refactor<double>(data_dir, refactor_dir, num_bitplanes, dims, ".dat", mode, qoi, approximator_eb, max_weight, block_size);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }

    clock_gettime(CLOCK_REALTIME, &end);
    double elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
    printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}
