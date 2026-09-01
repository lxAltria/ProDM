#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "ProDM/Reconstructor/MDR/Reconstructor.hpp"
#include "ProDM/Utils/QoIUtils.hpp"
#include "mpi.h"

using namespace std;

double local_reconstruct_time = 0;
unsigned long long local_retrieved_size = 0;
double local_max_abs_error = 0;
std::vector<uint32_t> offsets;
size_t metadata_size = 0;
double requested_tolerance = 0;
unsigned long long local_num_element = 0;

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
double compute_max_abs_error(const std::vector<T>& ori_data, const T * rec_data){
    double max_abs_error = 0;
    for(int i=0; i<ori_data.size(); i++){
        double abs_error = abs(ori_data[i] - rec_data[i]);
        if(abs_error > max_abs_error) max_abs_error = abs_error;
    }
    return max_abs_error;
}

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, vector<double>& tolerance, Reconstructor reconstructor){
    // struct timespec start, end;
    // int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    for(int i=0; i<tolerance.size(); i++){
        requested_tolerance = tolerance[i];
        // cout << "Start reconstruction" << endl;
        // err = clock_gettime(CLOCK_REALTIME, &start);
        local_reconstruct_time = -MPI_Wtime();
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i], -1);
        local_reconstruct_time += MPI_Wtime();
        // err = clock_gettime(CLOCK_REALTIME, &end);
        // cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        // auto dims = reconstructor.get_dimensions();
        local_max_abs_error = compute_max_abs_error(data, reconstructed_data);
        offsets = reconstructor.get_offsets();
        metadata_size = reconstructor.get_metadata_size();
        local_retrieved_size = static_cast<unsigned long long>(reconstructor.get_retrieved_size());
        // cout << "Retrieved data size = " << reconstructor.get_retrieved_size() << endl;
        // MGARD::print_statistics(data.data(), reconstructed_data, data.size(), retrieved_size);
        // std::cout << "Bitrate = " << (reconstructor.get_retrieved_size() * 8.0) / data.size() << std::endl;
        // COMP_UTILS::evaluate_gradients(data.data(), reconstructed_data, dims[0], dims[1], dims[2]);
        // COMP_UTILS::evaluate_average(data.data(), reconstructed_data, dims[0], dims[1], dims[2], 0);
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    local_num_element = static_cast<unsigned long long>(num_elements);
    // std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    T value_range = compute_global_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int argv_id = 1;
    string filename_base = string(argv[argv_id ++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string refactored_path_base = string(argv[argv_id++]);
    string output_path_base = string(argv[argv_id++]);

    int exp = static_cast<int>(std::round(std::log10(tolerance[0])));
	std::string wdata_file = output_path_base + "/1e" + std::to_string(exp) + ".bin";

    // int tmp_rank = rank;
    // int idx = tmp_rank / 64;
    // tmp_rank = tmp_rank % 64;
    // int idy = tmp_rank / 8;
    // tmp_rank = tmp_rank % 8;
    // int idz = tmp_rank;

    // string filename = filename_base + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + ".d64";
    // string refactored_path = refactored_path_base + "/" + to_string(idx) + "_" + to_string(idy) + "_" + to_string(idz) + "/";

    string filename = filename_base + to_string(rank) + ".d64";
    string refactored_path = refactored_path_base + "/" + to_string(rank) + "/";

    string metadata_file = refactored_path + "/metadata_mdr.bin";

    int num_levels = 0;
    int num_dims = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        // cout << "metadata_size = " << num_bytes << endl;
        // cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
    }
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = refactored_path + "/level_" + to_string(i) + "_mdr.bin";
        files.push_back(filename);
    }

    // using T = float;
    // using T_stream = uint32_t;
    using T = double;
    using T_stream = uint64_t;
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    auto encoder = MDR::PerBitBPEncoder_old<T, T_stream>();

    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();

    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    auto estimator = MDR::MaxErrorEstimatorHB<T>();
    auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    // auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);

    double global_reconstruct_time = 0;
    MPI_Reduce(&local_reconstruct_time, &global_reconstruct_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_elapsed_time = " << global_reconstruct_time << std::endl;

    if(!rank) std::cout << "Requested tolerance = " << requested_tolerance << std::endl;

    double global_max_abs_error = 0;
    MPI_Reduce(&local_max_abs_error, &global_max_abs_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_abs_error = " << global_max_abs_error << std::endl;

    unsigned long long global_retrieved_size = 0;
    MPI_Reduce(&local_retrieved_size, &global_retrieved_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

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
    size_t metadata_num_char = 0;
    auto metadata_data = MGARD::readfile<unsigned char>(metadata_file.c_str(), metadata_num_char);
    MPI_File_open(MPI_COMM_WORLD, wdata_file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &retrieved_file);
    MPI_File_write_at(retrieved_file, retrieved_offset, metadata_data.data(), metadata_data.size(), MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
    retrieved_offset += metadata_size;
    for(int i=0; i<offsets.size(); i++){
        size_t num_char = 0;
        auto level_data = MGARD::readfile<unsigned char>(files[i].c_str(), num_char);
        MPI_File_write_at(retrieved_file, retrieved_offset, level_data.data(), offsets[i], MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
        retrieved_offset += offsets[i];
    }
    MPI_File_close(&retrieved_file);
    MPI_Finalize();
    return 0;
}