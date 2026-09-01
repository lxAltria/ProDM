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

using namespace std;
bool negabinary = true;
bool greedy = false;
bool bfs = false;
bool greedy_bfs = false;
std::vector<int> coeff_interp_directions;
bool write_output = false;
std::string output_path;

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Start reconstruction, Requested Tolerance = " << tolerance[i] << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i], -1);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct_time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        auto dims = reconstructor.get_dimensions();
        size_t retrieved_size = reconstructor.get_retrieved_size();
        cout << "Retrieved data size = " << reconstructor.get_retrieved_size() << endl;
        MGARD::print_statistics(data.data(), reconstructed_data, data.size(), retrieved_size);
        std::cout << "Bitrate = " << (reconstructor.get_retrieved_size() * 8.0) / data.size() << std::endl;
        if(write_output) MGARD::writefile(output_path.c_str(), reconstructed_data, data.size());
        // COMP_UTILS::evaluate_gradients(data.data(), reconstructed_data, dims[0], dims[1], dims[2]);
        // COMP_UTILS::evaluate_average(data.data(), reconstructed_data, dims[0], dims[1], dims[2], 0);
    }
}

// Per level
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void init_fuse_composed_reconstructor(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::FuseComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = MDR::compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

// Per Layer
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void init_composed_reconstructor_new(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor_new<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = MDR::compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

// Per level + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void init_ordered_cp_reconstructor(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::OrderedCPReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    reconstructor.print();
    reconstructor.load_metadata();
    T value_range = MDR::compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

// // Per level + CP
// template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
// void init_ordered_cp_reconstructor(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
//     auto reconstructor = MDR::PartialOrderedCPReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
//     // reconstructor.print();
//     reconstructor.load_metadata();
//     T value_range = MDR::compute_value_range(data);
//     for(int i=0; i<tolerance.size(); i++){
//         tolerance[i] *= value_range;
//     }
//     evaluate(data, tolerance, reconstructor);
// }

// Per Layer + CP
template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void init_cp_reconstructor_new(std::vector<T>& data, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::CPReconstructor_new<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // reconstructor.print();
    reconstructor.load_metadata();
    T value_range = MDR::compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

void overall_reconstructing_initiator(std::string filename, std::string refactored_path, std::string data_type, std::vector<double>& tolerance, std::string interpreter_option, std::string encoder_option, std::string cp){
    std::string metadata_file = refactored_path + "/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    int interpolation_type = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
        if(!strcmp(data_type.c_str(), "-f")){
            interpolation_type = metadata[num_dims * sizeof(uint32_t) + num_levels * sizeof(float) + 2];
        } else if(!strcmp(data_type.c_str(), "-d")){
            interpolation_type = metadata[num_dims * sizeof(uint32_t) + num_levels * sizeof(double) + 2];
        } else {
            std::cerr << "Only two float type supported: -f or -d" << std::endl;
            exit(-1);
        }
    }
    std::cout << "interpolation_type = " << interpolation_type << std::endl;
    std::vector<std::string> files;
    for(int i=0; i<num_levels; i++){
        string filename = refactored_path + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    std::string data_file = refactored_path + "/data.bin";
    if(!strcmp(data_type.c_str(), "-f")){
        using T = float;
        using T_stream = uint32_t;
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        if (!interpolation_type){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
                // auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
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
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        }
    } else if(!strcmp(data_type.c_str(), "-d")){
        using T = double;
        using T_stream = uint64_t;
        size_t num_elements = 0;
        auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
        if (!interpolation_type){ // Per level
            auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
            auto interleaver = MDR::DirectInterleaver_new<T>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
            if(!strcmp(cp.c_str(), "-CP")){ // CP
                auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
                // auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_ordered_cp_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
                        init_fuse_composed_reconstructor(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
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
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_cp_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            } else { // no CP
                auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
                if(!strcmp(encoder_option.c_str(), "-Nega")){
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-XOR")){
                    auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                } else if (!strcmp(encoder_option.c_str(), "-PerBit")){
                    auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
                    if(!strcmp(interpreter_option.c_str(), "-BFS")){
                        auto interpreter = MDR::SignExcludeBFSBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-DP")){
                        auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    } else if (!strcmp(interpreter_option.c_str(), "-Greedy")){
                        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(estimator);
                        init_composed_reconstructor_new(data, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
                    }
                }
            }
        }
    } else {
        std::cerr << "Only two float type supported: -f or -d" << std::endl;
        exit(-1);
    }
}

void usage(char* cmd) {
    std::cout << "two_modes_reconstructor usage: " << cmd <<
                  " data_file -[dataType: f/d] num_of_tolerance tol1 tol2 ... toln refactored_path -[encoder_option: Nega/XOR/PerBit] -[interpreter_option: Greedy/DP/BFS] -[CP_or_not: CP/no_CP] [Optional: Reconstructed data path]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 -d 5 1e-1 1e-2 1e-3 1e-4 1e-5 /refactored/path -PerBit -BFS -CP" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string refactored_path = string(argv[argv_id++]);
    string data_type = string(argv[argv_id ++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    std::vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);
    }
    string encoder_option = string(argv[argv_id++]);
    string interpreter_option = string(argv[argv_id++]);
    string cp = string(argv[argv_id++]);
    if (argv_id < argc){
        output_path = string(argv[argv_id++]);
        write_output = true;
    }
    
    // std::cout << "argv_id = " << argv_id << ", argc = " << argc << std::endl;

    overall_reconstructing_initiator(filename, refactored_path, data_type, tolerance, interpreter_option, encoder_option, cp);
    
    return 0;
}