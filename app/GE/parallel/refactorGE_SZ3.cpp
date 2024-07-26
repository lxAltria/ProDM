#include <cmath>
#include <fstream>
#include <thread>
#include <chrono>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <zstd.h>
#include "mpi.h"
#include "adios2.h"
#include "SZ3/api/sz.hpp"
#include "recompose.hpp"

const std::vector<std::string> var_name{"U_aver", "V_aver", "W_aver", "Pressure", "Rho"};
const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
const int n_vars = 5;

template <class T>
void print_info(int p, const std::string& name, const std::vector<T>& vec){
    T max = vec[0];
    T min = vec[0];
    for(int i=1; i<vec.size(); i++){
        if(max < vec[i]) max = vec[i];
        if(min > vec[i]) min = vec[i];
    }
    printf("Processor %d var %s (size = %d): min = %.4f, max = %.4f\n", p, name.c_str(), vec.size(), min, max);
}

template <class T>
static T compute_vr(const std::vector<T>& vec){
    T min = vec[0];
    T max = vec[0];
    for(int i=0; i<vec.size(); i++){
        if(vec[i] < min) min = vec[i];
        if(vec[i] > max) max = vec[i];
    }
    return max - min;
}

template<class T>
char * SZ3_compress(size_t num_elements, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(num_elements);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
}

template<class Type>
void refactor_singleZone(const std::string data_prefix_path, int rank, int num_elements, std::vector<Type>& velocityX_vec, std::vector<Type>& velocityY_vec, std::vector<Type>& velocityZ_vec, std::vector<Type>& pressure_vec, std::vector<Type>& density_vec){
    std::string filename = data_prefix_path + "/data/sol_4114800_aver_b" + std::to_string(rank) + ".bp/";
    MGARD::writefile((filename + var_name_out[0] + ".dat").c_str(), velocityX_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[1] + ".dat").c_str(), velocityY_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[2] + ".dat").c_str(), velocityZ_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[3] + ".dat").c_str(), pressure_vec.data(), velocityX_vec.size());
    MGARD::writefile((filename + var_name_out[4] + ".dat").c_str(), density_vec.data(), velocityX_vec.size());
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){
            mask[i] = 1;
            num_valid_data ++;
        }
    }
    std::string mask_file = data_prefix_path + "/refactor/block_" + std::to_string(rank) + "_refactored/mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
    }
    std::vector<double> rel_ebs;
    const int num_snapshot = 10;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = data_prefix_path + "/refactor/block_" + std::to_string(rank) + "_refactored/" + var_name_out[i] + "/";
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_valid_data, buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
        }
        else{
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_elements, vars_vec[i].data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int argv_id = 1;
	std::string data_prefix_path(argv[argv_id++]);

	struct timespec start, end;
	int err;
	double elapsed_time;
	err = clock_gettime(CLOCK_REALTIME, &start);

    std::string filename = data_prefix_path + "/data/sol_4114800_aver_b" + std::to_string(rank) + ".bp";
    adios2::ADIOS ad; 
    adios2::IO reader_io = ad.DeclareIO("Input");
    adios2::Engine reader = reader_io.Open(filename, adios2::Mode::Read);

    while (true) {
        // Begin step
        adios2::StepStatus read_status = reader.BeginStep(adios2::StepMode::Read, 10.0f);
        if (read_status == adios2::StepStatus::NotReady) {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK) {
            break;
        }

        adios2::Variable<double> var_ad2;
        std::vector<std::vector<double>> vars_vec(n_vars, std::vector<double>());

        for (int i=0; i<n_vars; i++) {
            var_ad2 = reader_io.InquireVariable<double>(var_name[i]);
            std::vector<std::size_t> shape = var_ad2.Shape();
            reader.Get(var_ad2, vars_vec[i], adios2::Mode::Sync);
            reader.PerformGets();
            //print_info(rank, var_name[i], vars_vec[i]);
        }
	    refactor_singleZone(data_prefix_path, rank, vars_vec[0].size(), vars_vec[0], vars_vec[1], vars_vec[2], vars_vec[3], vars_vec[4]); 
        reader.EndStep();
    }
    reader.Close();

	err = clock_gettime(CLOCK_REALTIME, &end);
	printf("elapsed_time = %.6f\n", elapsed_time);

    MPI_Finalize();
    
    return 0;
}
