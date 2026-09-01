#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Reconstructor/MDR/Reconstructor.hpp"
#include "ProDM/Utils/QoIUtils.hpp"

using namespace std;

bool write_output = false;
string output_path;

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

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    T value_range = compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i]*value_range, -1);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
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

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    // tolerances are relative; evaluate() scales them by the value range
    evaluate(data, tolerance, reconstructor);
}

template <class T, class T_stream>
void launch_reconstructor(string filename, vector<double>& tolerance, int encoder_option, string metadata_file, const vector<string>& files){
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto encoder = MDR::XORNegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();

    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();

    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    auto estimator = MDR::MaxErrorEstimatorHB<T>();
    auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    // auto interpreter = MDR::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    if(encoder_option == 0){
        auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
    } else {
        auto encoder = MDR::PerBitBPEncoder_old<T, T_stream>();
        test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
    }
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_file refactored_dict num_tolerance tolerance1 ... toleranceN [Encoder: NegaBinary-0, PerBit-1] -[dataType: f/d] [Optional: reconstructed data path]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 3 1e-1 1e-2 1e-3 0 -d" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    string refactored_path = string(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    if(num_tolerance <= 0 || argc < argv_id + num_tolerance + 2){
        std::cerr << "Insufficient or invalid arguments (num_tolerance parsed as " << num_tolerance << "); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);
    }
    string metadata_file = refactored_path + "/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        if(num_bytes == 0){
            std::cerr << "Cannot read " << metadata_file << "; run the refactor first" << std::endl;
            return -1;
        }
        num_dims = metadata[0];
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        cout << "metadata_size = " << num_bytes << endl;
        cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
    }
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = refactored_path + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    if(argv_id + 1 >= argc){
        std::cerr << "Missing encoder option and/or data type option" << std::endl;
        usage(argv[0]);
        return -1;
    }
    int encoder_option = atoi(argv[argv_id ++]);
    std::string dtype = string(argv[argv_id++]);
    if(argv_id < argc){
        output_path = string(argv[argv_id ++]);
        write_output = true;
    }

    if (strcmp(dtype.c_str(), "-f") == 0){
        launch_reconstructor<float, uint32_t>(filename, tolerance, encoder_option, metadata_file, files);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        launch_reconstructor<double, uint64_t>(filename, tolerance, encoder_option, metadata_file, files);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }
    return 0;
}
