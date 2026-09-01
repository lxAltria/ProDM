#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "SZ3/api/sz.hpp"
#include "mpi.h"

using namespace std;

int file_ind = 0;
double local_elapsed_time = 0;

template<class T>
T compute_global_value_range(const std::vector<T>& data_vec){
	T global_max = 0, global_min = 0;
	T local_max = -std::numeric_limits<T>::max();
	T local_min = std::numeric_limits<T>::max();
	for(int i=0; i<data_vec.size(); i++){
		if(data_vec[i] > local_max) local_max = data_vec[i];
		if(data_vec[i] < local_min)	local_min = data_vec[i];
	}
	if(std::is_same<T, double>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	else if(std::is_same<T, float>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
	}
	return global_max - global_min;
}

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
}

inline int find_index(double target_rel_eb, double& rel_eb){
    int i = 0;
    while(target_rel_eb < rel_eb){
        i ++;
        rel_eb /= 10;
    }
    return i;
}

template<class T>
double compute_max_abs_error(const std::vector<T>& ori_data, const T * rec_data){
    double max_abs_error = 0;
    for(int i=0; i<ori_data.size(); i++){
        double abs_error = abs(ori_data[i] - rec_data[i]);
        if(abs_error > max_abs_error) max_abs_error = abs_error;
    }
    return max_abs_error;
}

template <class T>
T compute_value_range(const std::vector<T>& vec){
	T min = vec[0];
	T max = vec[0];
	for(int i=0; i<vec.size(); i++){
		if(vec[i] < min) min = vec[i];
		if(vec[i] > max) max = vec[i];
	}
	return max - min;
}

template<class T>
size_t PSZ3_delta_reconstructor(string refactor_dict, const vector<T>& data, double tolerance, std::vector<T>& reconstructed_data){
    T value_range = compute_value_range(data);
    double file_eb = 0.1;
    int current_ind = -1;
    file_ind = find_index(tolerance / value_range, file_eb);
    size_t retrieved_size = 0;
    T * tmp_reconstructed_data = (T *) malloc(data.size() * sizeof(T));
    if(file_ind > current_ind){
        for(int i=0; i<=file_ind; i++){
            string filename = refactor_dict + "/SZ3_delta_eb_" + std::to_string(i) + ".bin";
            size_t n = 0;
            auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
            retrieved_size += n;
            SZ3_decompress(cmpData.data(), n, tmp_reconstructed_data);
            int index = 0;
            for(int j=0; j<data.size(); j++){
                reconstructed_data[j] += tmp_reconstructed_data[j];
            }
        }
    }
    return retrieved_size;
}

template <class T>
void evaluate(string refactor_dict, string wdata_file, const vector<T>& data, std::vector<double>& tolerance){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // struct timespec start, end;
    // int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    T value_range = compute_global_value_range(data);
    size_t retrieved_size = 0;
    std::vector<T> reconstructed_data(data.size(), 0);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
        // cout << "Requested tolerance #" << i << " = " << tolerance[i] * value_range << endl;
        // cout << "Start reconstruction" << endl;
        // err = clock_gettime(CLOCK_REALTIME, &start);
        local_elapsed_time = -MPI_Wtime();
        retrieved_size = PSZ3_delta_reconstructor(refactor_dict, data, tolerance[i], reconstructed_data);
        local_elapsed_time += MPI_Wtime();
        // err = clock_gettime(CLOCK_REALTIME, &end);
        // cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        // cout << "Retrieved data size = " << retrieved_size << endl;
        // MGARD::print_statistics(data.data(), reconstructed_data.data(), data.size());
        // cout << "Bitrate = " << (retrieved_size * 8.0) / data.size() << std::endl;
        // cout << endl;
    }
    double global_elapsed_time = 0;
    MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_elapsed_time = " << global_elapsed_time << std::endl;

    if(!rank) std::cout << "Requested tolerance = " << tolerance[0] << std::endl;

    double local_max_abs_error = compute_max_abs_error(data, reconstructed_data.data());
    double global_max_abs_error = 0;
    MPI_Reduce(&local_max_abs_error, &global_max_abs_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_abs_error = " << global_max_abs_error << std::endl;

    unsigned long long local_retrieved_size = static_cast<unsigned long long>(retrieved_size);
    unsigned long long global_retrieved_size = 0;
    MPI_Reduce(&local_retrieved_size, &global_retrieved_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    unsigned long long local_num_element = static_cast<unsigned long long>(data.size());
    unsigned long long global_num_element = 0;
    MPI_Reduce(&local_num_element, &global_num_element, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if(!rank) std::cout << "Aggregated bitrate = " << global_retrieved_size * 8.0 / global_num_element << std::endl;

    unsigned long long retrieved_offset = 0;
    unsigned long long retrieved_buffer;

    for(int i=0; i<size; i++){
        if(i == rank){
            if(i != 0) {
                MPI_Recv(&retrieved_offset, 1, MPI_UNSIGNED_LONG_LONG, i-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            retrieved_buffer = retrieved_offset + local_retrieved_size;

            if(i != size - 1){
                MPI_Send(&retrieved_buffer, 1, MPI_UNSIGNED_LONG_LONG, i+1, 0, MPI_COMM_WORLD);
            }
        }
    }
    MPI_File retrieved_file;
    MPI_File_open(MPI_COMM_WORLD, wdata_file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &retrieved_file);
    for(int i=0; i<=file_ind; i++){
        string filename = refactor_dict + "/SZ3_delta_eb_" + std::to_string(i) + ".bin";
        size_t num_char = 0;
        auto level_data = MGARD::readfile<unsigned char>(filename.c_str(), num_char);
        MPI_File_write_at(retrieved_file, retrieved_offset, level_data.data(), num_char, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
        retrieved_offset += num_char;
    }
    MPI_File_close(&retrieved_file);
}

template <class T>
void test(string filename, string refactor_dict, string wdata_file, vector<double>& tolerance){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    // std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(refactor_dict, wdata_file, data, tolerance);
}

void usage(char* cmd) {
    std::cout << "PSZ3-delta usage: " << std::endl <<
                  "mpirun -n #num_cores para_PSZ3-delta_reconstructor data_file refactor_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d] output_path"
                  << std::endl
                  << "example: " << std::endl <<
                  "mpirun -n 2 para_PSZ3-delta_reconstructor density.d64 refactor/Density_refactored 3 1e-1 1e-2 1e-3 -d /path/output" << std::endl;
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if(!rank) usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string refactor_dict = string(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string dtype = string(argv[argv_id++]);
    string output_dict = string(argv[argv_id++]);

    // int tmp_rank = rank;
    // int idx = tmp_rank / 64;
    // tmp_rank = tmp_rank % 64;
    // int idy = tmp_rank / 8;
    // tmp_rank = tmp_rank % 8;
    // int idz = tmp_rank;

    // filename = filename + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + ".d64";
    // refactor_dict = refactor_dict + "/" + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + "/";

    filename = filename + to_string(rank) + ".d64";
    refactor_dict = refactor_dict + "/" + to_string(rank) + "/";

    int exp = static_cast<int>(std::round(std::log10(tolerance[0])));
	std::string wdata_file = output_dict + "/1e" + std::to_string(exp) + ".bin";

    if (strcmp(dtype.c_str(), "-f") == 0){
        using T = float;
        test<T>(filename, refactor_dict, wdata_file, tolerance);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        using T = double;
        test<T>(filename, refactor_dict, wdata_file, tolerance);
    }
    MPI_Finalize();
    return 0;
}