#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#define Dummy_Cmp 0
#define MGARD_Cmp 1
#define SZ2_Cmp 2
#define SZ3_Cmp 3
#define HPEZ_Cmp 4
#define GE_Cmp 5

using namespace std;

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
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    T value_range = compute_value_range(data);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Requested tolerance #" << i << " = " << tolerance[i] * value_range << endl;
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i] * value_range, -1);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        auto dims = reconstructor.get_dimensions();
        cout << "Retrieved data size = " << reconstructor.get_retrieved_size() << endl;
        MGARD::print_statistics(data.data(), reconstructed_data, data.size());
        cout << "Bitrate = " << (reconstructor.get_retrieved_size() * 8.0) / data.size() << std::endl;
        cout << endl;
    }
}

template <class T, class Approximator, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, const vector<double>& tolerance, Approximator approximator, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = PDR::ApproximationBasedReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(approximator, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(data, tolerance, reconstructor);
}

template <class T, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void launch_reconstructor(string filename, const vector<double>& tolerance, int approximator_rank, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    switch(approximator_rank){
        case Dummy_Cmp:{
            auto approximator = PDR::DummyApproximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        case MGARD_Cmp:{
            auto approximator = PDR::MGARDApproximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        case SZ2_Cmp:{
            auto approximator = PDR::SZ2Approximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        case SZ3_Cmp:{
            auto approximator = PDR::SZ3Approximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        case HPEZ_Cmp:{
            auto approximator = PDR::HPEZApproximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        case GE_Cmp:{
            auto approximator = PDR::GEApproximator<T>();
            test<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
            break;
        }
        default:
            perror("Undefined Approximator\n");
            break;
    }
}

void usage(char* cmd) {
    std::cout << "QPro usage: " << cmd <<
                  " data_file refactored_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactored/Density_refactored 3 1e-1 1e-2 1e-3 -d 4" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    std::string filename = string(argv[argv_id++]);
    std::string refactor_dict = string(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    string metadata_file = refactor_dict + "/metadata.bin";
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
        string filename = refactor_dict + "/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

    using T_stream = uint32_t;

    std::string dtype = string(argv[argv_id++]);
    int approximator = atoi(argv[argv_id++]);
    if (strcmp(dtype.c_str(), "-f") == 0){
        using T = float;
        auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
        auto estimator = MDR::MaxErrorEstimatorHB<T>();
        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
        launch_reconstructor<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        using T = double;
        auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
        auto estimator = MDR::MaxErrorEstimatorHB<T>();
        auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
        launch_reconstructor<T>(filename, tolerance, approximator, encoder, compressor, estimator, interpreter, retriever);
    }

    return 0;
}