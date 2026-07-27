#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "nomask_Synthesizer4GE.hpp"

using namespace MDR;

int main(int argc, char** argv){

    using T = float;
    int argv_id = 1;
    std::string data = argv[argv_id++];
    std::cout << "data = " << data << std::endl;
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";
	int max_weight = 4;
	int block_size = 1;
	if(argc > 3){
		max_weight = atoi(argv[argv_id++]);
		block_size = atoi(argv[argv_id++]);
	}

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    if(std::strcmp(data.c_str(), "Hurricane") == 0){
		uint32_t dims[3] = {100, 500, 500};
        refactor_velocities_3D_PMGARD_WBP<T>(data, dims[0], dims[1], dims[2], data_file_prefix, rdata_file_prefix, max_weight, block_size);
    }
	else if (std::strcmp(data.c_str(), "NYX") == 0){
		refactor_velocities_3D_PMGARD_WBP<T>(data, 512, 512, 512, data_file_prefix, rdata_file_prefix, max_weight, block_size);
	}
	else if (std::strcmp(data.c_str(), "GE_small") == 0){
		refactor_velocities_1D_PMGARD_WBP<T>(data_file_prefix, rdata_file_prefix, max_weight, block_size);
	}
    
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;

}