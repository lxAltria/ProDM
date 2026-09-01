#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#include "PDR/Reconstructor/ApproximationBasedDeltaReconstructor.hpp"
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

template <class T, class Approximator>
void test(string filename, string refactor_dict, const vector<double>& tolerance, Approximator approximator){
    auto reconstructor = PDR::ApproximationBasedDeltaReconstructor<T, Approximator>(approximator, refactor_dict);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(data, tolerance, reconstructor);
}

template <class T>
void launch_reconstructor(string filename, string refactor_dict, const vector<double>& tolerance, int approximator_rank){
    switch(approximator_rank){
        case Dummy_Cmp:{
            auto approximator = PDR::DummyApproximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        case MGARD_Cmp:{
            auto approximator = PDR::MGARDApproximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        case SZ2_Cmp:{
            auto approximator = PDR::SZ2Approximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        case SZ3_Cmp:{
            auto approximator = PDR::SZ3Approximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        case HPEZ_Cmp:{
            auto approximator = PDR::HPEZApproximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        case GE_Cmp:{
            auto approximator = PDR::GEApproximator<T>();
            test<T>(filename, refactor_dict, tolerance, approximator);
            break;
        }
        default:
            perror("Undefined Approximator\n");
            break;
    }
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_file refactored_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4, GE-5]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 3 1e-1 1e-2 1e-3 -d 4" << std::endl;
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

    if(argv_id + 1 >= argc){
        std::cerr << "Missing data type option and/or approximator option" << std::endl;
        usage(argv[0]);
        return -1;
    }
    std::string dtype = string(argv[argv_id++]);
    int approximator = atoi(argv[argv_id++]);

    {
        std::string metadata_file = refactor_dict + "/metadata.bin";
        FILE * file = fopen(metadata_file.c_str(), "r");
        if(file == NULL){
            std::cerr << "Cannot open " << metadata_file << "; run test_pdr_delta_refactor first" << std::endl;
            return -1;
        }
        fclose(file);
    }

    if (strcmp(dtype.c_str(), "-f") == 0){
        launch_reconstructor<float>(filename, refactor_dict, tolerance, approximator);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        launch_reconstructor<double>(filename, refactor_dict, tolerance, approximator);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }

    return 0;
}
