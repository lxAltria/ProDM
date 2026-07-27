#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "qoi_utils.hpp"

using namespace std;

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i], -1);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        auto dims = reconstructor.get_dimensions();
        size_t retrieved_size = reconstructor.get_retrieved_size();
        cout << "Retrieval size = " << reconstructor.get_retrieved_size() << endl;
        MGARD::print_statistics(data.data(), reconstructed_data, data.size(), retrieved_size);
        std::cout << "Bitrate = " << (reconstructor.get_retrieved_size() * 8.0) / data.size() << std::endl;
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::OrderedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    T value_range = MDR::compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        tolerance[i] *= value_range;
    }
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string data_type = string(argv[argv_id ++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string refactored_path = string(argv[argv_id ++]);
    int encoder_option = atoi(argv[argv_id ++]);
    string metadata_file = refactored_path + "/refactored_data/metadata.bin";
    string data_file = refactored_path + "/refactored_data/data.bin";
    int num_levels = 0;
    int num_dims = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
    }

    if(!strcmp(data_type.c_str(), "-f")){
        using T = float;
        using T_stream = uint32_t;
        // using T = double;
        // using T_stream = uint64_t;
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        // auto decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
        auto interleaver = MDR::DirectInterleaver<T>();
        // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
        auto estimator = MDR::MaxErrorEstimatorHB<T>();
        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
        if(encoder_option == 0){
            auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        } else if (encoder_option == 1) {
            auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }else {
            auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }
    } else if (!strcmp(data_type.c_str(), "-d")){
        using T = double;
        using T_stream = uint64_t;
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        // auto decomposer = MDR::MGARDHierarchical_Coeff_Decomposer_Interleaver<T>(0);
        auto interleaver = MDR::DirectInterleaver<T>();
        // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto retriever = MDR::OrderedFileRetriever(metadata_file, data_file);
        auto estimator = MDR::MaxErrorEstimatorHB<T>();
        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
        if(encoder_option == 0){
            auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        } else if (encoder_option == 1) {
            auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }else {
            auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }
    } else {
        std::cerr << "Only two float type supported: -f or -d" << std::endl;
        exit(-1);
    }
    return 0;
}