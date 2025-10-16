#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "nomask_Synthesizer4GE.hpp"
#define Dummy_Cmp 0
#define SZ3_Cmp 1
#define PMGARD 2
#define GE_Cmp 3
#define HPEZ_Cmp 4
#define SZ2_Cmp 5
#define MGARD_Cmp 6
#define HPEZ_Vtot2 7
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
    T approximator_eb = 0.001;
	if(argc > 5){
		max_weight_for_vtot = atoi(argv[argv_id++]);
        if (std::strcmp(data.c_str(), "GE") == 0){
            max_weight_for_temperature = atoi(argv[argv_id++]);
        }
        else{
		    block_size = atoi(argv[argv_id++]);
        }
        approximator_eb = atof(argv[argv_id++]);
	}

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    switch (compressor)
    {
    case Dummy_Cmp:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_Dummy_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_Dummy_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_Dummy_BP<T>(data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_Dummy_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_Dummy_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_Dummy_WBP<T>(data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, max_weight_for_temperature, block_size);
            }
        }
        break;
    case SZ3_Cmp:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 100*500*500);
                refactor_velocities_3D_SZ3_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 512*512*512);
                refactor_velocities_3D_SZ3_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 98*1200*1200);
                refactor_velocities_3D_SZ3_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                // refactor_velocities_A1D_SZ3_BP<T>(data_file_prefix, rdata_file_prefix, 256*384*384);
                refactor_velocities_3D_SZ3_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ3_BP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_SZ3_BP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "GE") == 0){
                refactor_velocities_1D_SZ3_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight_for_vtot, approximator_eb, max_weight_for_temperature, block_size);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_SZ3_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
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
    case GE_Cmp:
        if(!weighted){
            refactor_velocities_1D_GE_BP<T>(data_file_prefix, rdata_file_prefix, approximator_eb);
        }
        else{
            refactor_velocities_1D_GE_WBP<T>(data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, max_weight_for_temperature);
        }
        break;
    case HPEZ_Cmp:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_S3D_xixj_BP<T>(data, 1200, 334, 200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_S3D_xixj_WBP<T>(data, 1200, 334, 200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_HPEZ_WBP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
        }
        break;
    case SZ2_Cmp:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_SZ2_BP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_SZ2_WBP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
        }
        break;
    case MGARD_Cmp:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_MGARD_BP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 500, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_MGARD_WBP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
            }
        }
        break;
    case HPEZ_Vtot2:
        if(!weighted){
            if(std::strcmp(data.c_str(), "Hurricane") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "NYX") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "SCALE") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Miranda") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "S3D") == 0){
                refactor_S3D_xixj_BP<T>(data, 1200, 334, 200, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
            else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                refactor_velocities_3D_HPEZ_BP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb);
            }
        }
        else{
            int max_aggregated_weight = atoi(argv[argv_id++]);
            if(!max_aggregated_weight){
                if(std::strcmp(data.c_str(), "Hurricane") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "NYX") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "SCALE") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "Miranda") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "S3D") == 0){
                    refactor_S3D_xixj_WBP<T>(data, 1200, 334, 200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                    refactor_velocities_square_3D_HPEZ_WBP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
            }
            else{
                if(std::strcmp(data.c_str(), "Hurricane") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 100, 500, 500, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "NYX") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "SCALE") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 98, 1200, 1200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "Miranda") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 256, 384, 384, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "S3D") == 0){
                    refactor_S3D_xixj_WBP<T>(data, 1200, 334, 200, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "Nek5000") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 510, 510, 510, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "JHTDB_3GB") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
                else if (std::strcmp(data.c_str(), "JHTDB_1.5GB") == 0){
                    refactor_velocities_and_square_3D_HPEZ_WBP<T>(data, 256, 512, 512, data_file_prefix, rdata_file_prefix, approximator_eb, max_weight_for_vtot, block_size);
                }
            }
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