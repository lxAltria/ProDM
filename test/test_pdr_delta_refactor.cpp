#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <sys/stat.h>
#include "ProDM/MGARDx/utils.hpp"
#include "ProDM/Refactor/PDR/Refactor.hpp"
#include "ProDM/Refactor/PDR/ApproximationBasedDeltaRefactor.hpp"
#define Dummy_Cmp 0
#define MGARD_Cmp 1
#define SZ2_Cmp 2
#define SZ3_Cmp 3
#define HPEZ_Cmp 4
#define GE_Cmp 5

using namespace std;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    // target_level and num_bitplanes are not used by the delta refactor
    refactor.refactor(data.data(), dims, 0, 0);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Approximator>
void test(string filename, string refactor_dict, const vector<uint32_t>& dims, Approximator approximator){
    auto refactor = PDR::ApproximationBasedDeltaRefactor<T, Approximator>(approximator, refactor_dict);
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, refactor);
}

template <class T>
void launch_refactor(string filename, string refactor_dict, const vector<uint32_t>& dims, int approximator_rank){
    switch(approximator_rank){
        case Dummy_Cmp:{
            auto approximator = PDR::DummyApproximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#ifdef PRODM_HAVE_MGARD
        case MGARD_Cmp:{
            auto approximator = PDR::MGARDApproximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#endif
#ifdef PRODM_HAVE_SZ2
        case SZ2_Cmp:{
            auto approximator = PDR::SZ2Approximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#endif
#ifdef PRODM_HAVE_SZ3
        case SZ3_Cmp:{
            auto approximator = PDR::SZ3Approximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#endif
#ifdef PRODM_HAVE_HPEZ
        case HPEZ_Cmp:{
            auto approximator = PDR::HPEZApproximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#endif
#ifdef PRODM_HAVE_HPEZ
        case GE_Cmp:{
            auto approximator = PDR::GEApproximator<T>();
            test<T>(filename, refactor_dict, dims, approximator);
            break;
        }
#endif
        default:
            std::cerr << "Approximator " << approximator_rank << " is unknown or not enabled at build time (see PRODM_WITH_* CMake options)" << std::endl;
            break;
    }
}

void usage(char* cmd) {
    std::cout << "usage: " << cmd <<
                  " data_file output_dict num_dim dim0 .. dimn -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4, GE-5]"
                  << std::endl
                  << "example: " << cmd <<
                  " density.d64 refactor/Density_refactored 3 256 384 384 -d 4" << std::endl;
}

int main(int argc, char ** argv){
    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }
    int argv_id = 1;
    std::string filename = string(argv[argv_id ++]);
    std::string refactor_dict = string(argv[argv_id ++]);
    int num_dims = atoi(argv[argv_id ++]);
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id ++]);
    }

    if(argv_id + 1 >= argc){
        std::cerr << "Missing data type option and/or approximator option" << std::endl;
        usage(argv[0]);
        return -1;
    }
    std::string dtype = string(argv[argv_id++]);
    int approximator = atoi(argv[argv_id++]);

    mkdir(refactor_dict.c_str(), 0777);

    if (strcmp(dtype.c_str(), "-f") == 0){
        launch_refactor<float>(filename, refactor_dict, dims, approximator);
    } else if (strcmp(dtype.c_str(), "-d") == 0){
        launch_refactor<double>(filename, refactor_dict, dims, approximator);
    } else {
        std::cerr << "Unknown data type option: " << dtype << " (expected -f or -d); check the argument order" << std::endl;
        usage(argv[0]);
        return -1;
    }

    return 0;
}
