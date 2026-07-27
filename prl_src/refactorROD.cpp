#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <mpi.h>
#include "nomask_Synthesizer4GE.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "MDR/Refactor/Refactor.hpp"
#include "SZ3/api/sz.hpp"
#include "PDR/Refactor/Refactor.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#include "ompSZp_typemanager.h"

int max_weight_for_vtot = 4;
int max_weight_for_pressure = 3;
int max_weight_for_density = 3;
const std::string data_prefix_path = "/pscratch/xli281_uksr/xliang/GE";

template<class T>
void refactor_singleZone(int rank){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityX.dat.rod").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityY.dat.rod").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityZ.dat.rod").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/Pressure.dat.rod").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/Density.dat.rod").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    std::vector<int> max_weights = {max_weight_for_vtot, max_weight_for_vtot, max_weight_for_vtot, max_weight_for_pressure, max_weight_for_density};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);

    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    std::vector<std::vector<T>> Vtot(3, std::vector<T>(num_elements));
    std::vector<std::vector<T>> Temp(2, std::vector<T>(num_elements));
    T R = 287.1;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data++;
        }
        if(mask[i]){
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        Temp[0][i] = 1.0 / pow(density_vec[i], 2.0);
        Temp[1][i] = 1.0 / pow(density_vec[i], 2.0);
    }

    std::string mask_file = data_prefix_path + "/block_" + std::to_string(rank) + "_refactored/mask.bin";
    // std::string mask_file = data_prefix_path + "/additional/block_" + std::to_string(rank) + "_mask.bin";
    writemask(mask_file.c_str(), mask.data(), mask.size());
    uint8_t * weight_data = (uint8_t *) malloc(num_elements);
    std::string block_path = data_prefix_path + "/block_" + std::to_string(rank) + "_refactored/block_sizes.dat";
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = data_prefix_path + "/block_" + std::to_string(rank) + "_refactored/" + var_list[i] + "/";
        std::string metadata_file = rdir_prefix + "metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::GEApproximator<T>();
        auto compressor = MDR::AdaptiveLevelCompressor(32);
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint64_t>();
        bool negabinary = true;
        // auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        // bool negabinary = true;
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateWBPRefactor_GE<T>(approximator, encoder, compressor, writer, negabinary);
        if(i < 3){
            refactor.QoI = Vtot[i];
            refactor.mask = mask;
        }
        else refactor.QoI = Temp[i-3];

        std::string approximator_path = data_prefix_path + "/block_" + std::to_string(rank) + "_refactored/" + var_list[i] + "/approximator.dat";
        // std::string approximator_path = data_prefix_path + "/additional/block_" + std::to_string(rank) + "_" + var_list[i] + "_approximator.dat";
        refactor.refactor(vars_vec[i].data(), weight_data, block_path, approximator_path, dims, target_level, num_bitplanes, max_weights[i]);  
        std::string weight_file = data_prefix_path + "/block_" + std::to_string(rank) + "_refactored/" + var_list[i] + "/weight.bin";
        // std::string weight_file = data_prefix_path + "/additional/block_" + std::to_string(rank) + "_" + var_list[i] + "_weight.bin";
        MGARD::writefile(weight_file.c_str(), weight_data, refactor.get_weight_size());
        // std::cout << "rank = " << rank << ", i = " << i << ", get_weight_size = " << refactor.get_weight_size() << std::endl;
    }
    free(weight_data);
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);
    refactor_singleZone<double>(rank);
	err = clock_gettime(CLOCK_REALTIME, &end);
    elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

    double max_time;
	MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(rank==0) printf("elapsed_time = %.6f\n", max_time);

    MPI_Finalize();
    
    return 0;
}