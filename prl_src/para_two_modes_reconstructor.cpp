#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "ProDM/MGARDx/utils.hpp"
#include "ProDM/Reconstructor/MDR/Reconstructor.hpp"
#include "ProDM/Utils/QoIUtils.hpp"
#include "mpi.h"

using namespace std;
bool negabinary = true;
bool greedy = false;
bool bfs = false;
bool greedy_bfs = false;
std::vector<int> coeff_interp_directions;
std::vector<uint32_t> offsets;
double local_max_abs_error = 0;
std::string metadata_file;
std::vector<std::string> files;
std::string data_file;
uint32_t metadata_size = 0;
unsigned long long local_retrieved_size = 0;
unsigned long long local_num_element = 0;
double requested_tolerance = 0;

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
double evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    double total_reconstruct_time = 0;
    for(int i=0; i<tolerance.size(); i++){
        requested_tolerance = tolerance[i];
        double start = MPI_Wtime();
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i], -1);
        total_reconstruct_time = MPI_Wtime() - start;
        local_max_abs_error = compute_max_abs_error(data, reconstructed_data);
        offsets = reconstructor.get_offsets();
        metadata_size = reconstructor.get_metadata_size();
        local_retrieved_size = static_cast<unsigned long long>(reconstructor.get_retrieved_size());
        // auto dims = reconstructor.get_dimensions();
        // size_t retrieved_size = reconstructor.get_retrieved_size();
        // MGARD::print_statistics(data.data(), reconstructed_data, data.size(), retrieved_size);
        // std::cout << "Bitrate = " << (reconstructor.get_retrieved_size() * 8.0) / data.size() << std::endl;
    }
    return total_reconstruct_time;
}

// Per level
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
double init_fuse_composed_reconstructor(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::FuseComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = compute_global_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    return evaluate(data, tolerance, reconstructor);
}

// Per Layer
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
double init_composed_reconstructor_new(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor_new<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = compute_global_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    return evaluate(data, tolerance, reconstructor);
}

// Per level + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
double init_ordered_cp_reconstructor(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::OrderedCPReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = compute_global_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    return evaluate(data, tolerance, reconstructor);
}

// Per Layer + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
double init_cp_reconstructor_new(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::CPReconstructor_new<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = compute_global_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    return evaluate(data, tolerance, reconstructor);
}

double overall_reconstructing_initiator(std::string filename, std::string refactored_path, std::string data_type, std::vector<double>& tolerance, std::string interpreter_option, std::string encoder_option, std::string cp){
    metadata_file = refactored_path + "/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    int interpolation_type = 0;
    double reconstruct_time = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        // cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
        if(!strcmp(data_type.c_str(), "-f")){
            interpolation_type = metadata[num_dims * sizeof(uint32_t) + num_levels * sizeof(float) + 2];
        } else if(!strcmp(data_type.c_str(), "-d")){
            interpolation_type = metadata[num_dims * sizeof(uint32_t) + num_levels * sizeof(double) + 2];
        } else {
            std::cerr << "Only two float type supported: -f or -d" << std::endl;
            exit(-1);
        }
    }
    // std::cout << "interpolation_type = " << interpolation_type << std::endl;
    for(int i=0; i<num_levels; i++){
        string filename = refactored_path + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    data_file = refactored_path + "/data.bin";
    if(!strcmp(data_type.c_str(), "-f")){
        using T = float;
        using T_stream = uint32_t;
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        local_num_element = static_cast<unsigned long long>(num_elements);
        if (!interpolation_type){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        } else { // Per Layer
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        }
    } else if(!strcmp(data_type.c_str(), "-d")){
        using T = double;
        using T_stream = uint64_t;
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        local_num_element = static_cast<unsigned long long>(num_elements);
        if (!interpolation_type){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        reconstruct_time = init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        } else { // Per Layer
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        reconstruct_time = init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        }
    } else {
        std::cerr << "Only two float type supported: -f or -d" << std::endl;
        exit(-1);
    }
    return reconstruct_time;
}

void usage(char* cmd) {
    std::cout << "two_modes_reconstructor usage: " << cmd <<
                  " data_file -[dataType: f/d] num_of_tolerance tol1 tol2 ... toln refactored_path -[encoder_option: Nega/XOR/PerBit] -[interpreter_option: Greedy/DP/BFS] -[CP_or_not: CP/no_CP] output_path"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 -d 5 1e-1 1e-2 1e-3 1e-4 1e-5 /refactored/path -PerBit -eb -CP /output/path" << std::endl;
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if(rank == 0) usage(argv[0]);
        MPI_Finalize();
        return 0;
    }
    int argv_id = 1;
    string filename_base = string(argv[argv_id ++]);
    string data_type = string(argv[argv_id ++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    std::vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);
    }
    string refactored_path_base = string(argv[argv_id++]);
    string encoder_option = string(argv[argv_id++]);
    string interpreter_option = string(argv[argv_id++]);
    string cp = string(argv[argv_id++]);
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

    double local_reconstruct_time = overall_reconstructing_initiator(filename, refactored_path, data_type, tolerance, interpreter_option, encoder_option, cp);
    // cout << "[Rank " << rank << "] Reconstruct time: " << local_reconstruct_time << "s" << endl;

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
    if(offsets.size() == 1){
        size_t num_char = 0;
        auto level_data = MGARD::readfile<unsigned char>(data_file.c_str(), num_char);
        MPI_File_write_at(retrieved_file, retrieved_offset, level_data.data(), offsets[0], MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
    } else {
        for(int i=0; i<offsets.size(); i++){
            size_t num_char = 0;
            auto level_data = MGARD::readfile<unsigned char>(files[i].c_str(), num_char);
            MPI_File_write_at(retrieved_file, retrieved_offset, level_data.data(), offsets[i], MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
            retrieved_offset += offsets[i];
        }
    }
    MPI_File_close(&retrieved_file);
    MPI_Finalize();
    return 0;
}