#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"

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
        cout << "Retrieved data size = " << reconstructor.get_retrieved_size() << endl;
        MGARD::print_statistics(data.data(), reconstructed_data, data.size(), retrieved_size);
        // COMP_UTILS::evaluate_gradients(data.data(), reconstructed_data, dims[0], dims[1], dims[2]);
        // COMP_UTILS::evaluate_average(data.data(), reconstructed_data, dims[0], dims[1], dims[2], 0);
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, const vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::FuseComposedReconstructor_new<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string refactored_path = string(argv[argv_id++]);
    string metadata_file = refactored_path + "/refactored_data/metadata.bin";
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
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = refactored_path + "/refactored_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

    // using T = float;
    // using T_stream = uint32_t;
    using T = double;
    using T_stream = uint64_t;
    auto decomposer = MDR::MGARDHierarchicalDecomposer_Interleaver_new<T>();
    auto interleaver = MDR::DirectInterleaver_new<T>();
    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();

    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();

    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    auto estimator = MDR::MaxErrorEstimatorHB<T>();
    auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    // auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    // auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
    // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(estimator);
    test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
    return 0;
}