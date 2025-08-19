#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "nomask_Synthesizer4GE.hpp"
#define Dummy 0
#define SZ3 1
#define PMGARD 2
#define GE 3
using namespace MDR;

int main(int argc, char** argv){

    using T = double;
    int argv_id = 1;
    int compressor = atoi(argv[argv_id++]);
    int weighted = atoi(argv[argv_id++]);
    std::string data = argv[argv_id++];
    // std::cout << "data = " << data << std::endl;
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";
    int max_weight_for_vtot = 4;
    int max_weight_for_temperature = 0;
	int block_size = 1;
	if(argc > 5){
		max_weight_for_vtot = atoi(argv[argv_id++]);
        if (std::strcmp(data.c_str(), "GE") == 0){
            max_weight_for_temperature = atoi(argv[argv_id++]);
        }
        else{
		    block_size = atoi(argv[argv_id++]);
        }
	}

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    switch (compressor)
    {
    case Dummy:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_Dummy_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_Dummy_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_Dummy_BP<T>(data_file_prefix, rdata_file_prefix);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_Dummy_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_Dummy_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_Dummy_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight_for_vtot, max_weight_for_temperature, block_size);
            }
        }
        break;
    case SZ3:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 100*500*500);
                refactor_velocities_3D_SZ3_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 512*512*512);
                refactor_velocities_3D_SZ3_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 98*1200*1200);
                refactor_velocities_3D_SZ3_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 256*384*384);
                refactor_velocities_3D_SZ3_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_SZ3_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight_for_vtot, max_weight_for_temperature, block_size);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ3_WBP_JHTDB<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
        }
        break;
    case PMGARD:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                // refactor_velocities_A1D_PMGARD_BP<T>(data_file_prefix, rdata_file_prefix, 100*500*500);
                refactor_velocities_3D_PMGARD_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                // refactor_velocities_A1D_PMGARD_BP<T>(data_file_prefix, rdata_file_prefix, 512*512*512);
                refactor_velocities_3D_PMGARD_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_PMGARD_BP<T>(data_file_prefix, rdata_file_prefix);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_PMGARD_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_PMGARD_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_PMGARD_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight_for_vtot, max_weight_for_temperature, block_size);
            }
        }
        break;
    case GE:
        if(!weighted){
            refactor_velocities_1D_GE_BP<T>(data_file_prefix, rdata_file_prefix);
        }
        else{
            refactor_velocities_1D_GE_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight_for_vtot, max_weight_for_temperature);
        }
        break;
    default:
        break;
    }
    
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;

}