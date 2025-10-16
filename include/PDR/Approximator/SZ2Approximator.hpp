#ifndef _PDR_SZ2_APPROXIMATOR_HPP
#define _PDR_SZ2_APPROXIMATOR_HPP

#include "ApproximatorInterface.hpp"
#include "sz/sz_api.h"
#include "utils.hpp"

namespace PDR {
    // SZ approximator with SZ2 prediction
    template<class T>
    class SZ2Approximator : public concepts::ApproximatorInterface<T> {
    public:
        SZ2Approximator(){}
        size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximator_eb, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            if (SZ_NSCS == SZ_Init(NULL)) {
                printf("SZ initialization failed!\n");
                exit(0);
            }
            size_t cmpSize = 0;
            unsigned char *cmpData = NULL;
            size_t num_elements = std::accumulate(dimensions.begin(), dimensions.end(), size_t{1}, std::multiplies<size_t>());
            int datatype = std::is_same<T, double>::value ? SZ_DOUBLE : SZ_FLOAT;
            confparams_cpr->errorBoundMode = ABS;
            confparams_cpr->absErrBound = approximator_eb;
            size_t r5 = 0;
            size_t r4 = 0;
            size_t r3 = 0;
            size_t r2 = 0;
            size_t r1 = 0;
            if(dimensions.size() == 1){
                r1 = dimensions[0];
            }
            else if(dimensions.size() == 2){
                r1 = dimensions[0];
                r2 = dimensions[1];
            }
            else if(dimensions.size() == 3){
                r1 = dimensions[0];
                r2 = dimensions[1];
                r3 = dimensions[2];
            }
            else{
                std::cout << "Dimension larger than 4 is not supported at SZ approximator!" << std::endl;
                exit(-1);
            }
            cmpData = ::SZ_compress(datatype, data, &cmpSize, r5, r4, r3, r2, r1);
            if(filename.size()) approximator_file_name = filename;
            MGARD::writefile(approximator_file_name.c_str(), cmpData, cmpSize);
            approximator_file_size = cmpSize;
            // std::cout << "Approximator size = " << approximator_file_size << std::endl;
            // std::cout << "num_elements = " << num_elements << std::endl;
            T * dec_data = (T *) malloc(num_elements * sizeof(T));
            dec_data = static_cast<T*>(::SZ_decompress(datatype, cmpData, cmpSize, r5, r4, r3, r2, r1));
            free(cmpData);
            // T tmp = 0;
            for(int i=0; i<num_elements; i++){
                data[i] -= dec_data[i];
            //     if(fabs(data[i]) > tmp) tmp = fabs(data[i]);
            }
            // std::cout << "max diff = " << tmp << std::endl;
            free(dec_data);
            SZ_Finalize();
            return approximator_file_size;
        }

        void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            if(filename.size()) approximator_file_name = filename;
            size_t num = 0;
            auto cmpData = MGARD::readfile<unsigned char>(approximator_file_name.c_str(), num);
            approximator_file_size = num;
            int datatype = std::is_same<T, double>::value ? SZ_DOUBLE : SZ_FLOAT;
            size_t r5 = 0;
            size_t r4 = 0;
            size_t r3 = 0;
            size_t r2 = 0;
            size_t r1 = 0;
            if(dimensions.size() == 1){
                r1 = dimensions[0];
            }
            else if(dimensions.size() == 2){
                r1 = dimensions[0];
                r2 = dimensions[1];
            }
            else if(dimensions.size() == 3){
                r1 = dimensions[0];
                r2 = dimensions[1];
                r3 = dimensions[2];
            }
            else{
                std::cout << "Dimension larger than 4 is not supported at SZ approximator!" << std::endl;
                exit(-1);
            }
            T* decData = static_cast<T*>(::SZ_decompress(datatype, cmpData.data(), num, r5, r4, r3, r2, r1));
            size_t num_elements = std::accumulate(dimensions.begin(), dimensions.end(), size_t{1}, std::multiplies<size_t>());
            memcpy(data, decData, num_elements * sizeof(T));
            free(decData);
        }

        size_t get_size() const {
            return approximator_file_size;
        }

        void print() const {
            std::cout << "PDR SZ2 approximator" << std::endl;
        }
    private:
        size_t approximator_file_size;
        std::string approximator_file_name = "approximator.dat";
    };
}

#endif