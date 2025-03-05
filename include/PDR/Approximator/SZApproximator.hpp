#ifndef _PDR_SZ_APPROXIMATOR_HPP
#define _PDR_SZ_APPROXIMATOR_HPP

#include "ApproximatorInterface.hpp"
#include "SZ3/api/sz.hpp"
#include "utils.hpp"

namespace PDR {
    // SZ approximator with SZ3 prediction
    template<class T>
    class SZApproximator : public concepts::ApproximatorInterface<T> {
    public:
        SZApproximator(){}
        size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximator_eb, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {

            size_t cmpSize = 0;
            char *cmpData = NULL;
            size_t num_elements = 0;
            if(dimensions.size() == 1){
                SZ3::Config conf(dimensions[0]);
                conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
                conf.errorBoundMode = SZ3::EB_ABS;
                conf.absErrorBound = approximator_eb;
                cmpData = SZ_compress<T>(conf, data, cmpSize);
                num_elements = conf.num;
            }
            else if(dimensions.size() == 2){
                SZ3::Config conf(dimensions[0], dimensions[1]);
                conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
                conf.errorBoundMode = SZ3::EB_ABS;
                conf.absErrorBound = approximator_eb;
                cmpData = SZ_compress<T>(conf, data, cmpSize);
                num_elements = conf.num;
            }
            else if(dimensions.size() == 3){
                SZ3::Config conf(dimensions[0], dimensions[1], dimensions[2]);
                conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
                conf.errorBoundMode = SZ3::EB_ABS;
                conf.absErrorBound = approximator_eb;
                cmpData = SZ_compress<T>(conf, data, cmpSize);
                num_elements = conf.num;
            }
            else{
                std::cout << "Dimension larger than 4 is not supported at SZ approximator!" << std::endl;
                exit(-1);
            }
            if(filename.size()) approximator_file_name = filename;
            MGARD::writefile(approximator_file_name.c_str(), cmpData, cmpSize);
            approximator_file_size = cmpSize;
            // std::cout << "Approximator size = " << approximator_file_size << std::endl;
            // std::cout << "num_elements = " << num_elements << std::endl;
            T * dec_data = (T *) malloc(num_elements * sizeof(T));
            SZ3::Config conf;
            SZ_decompress<T>(conf, cmpData, cmpSize, dec_data);
            free(cmpData);
            T tmp = 0;
            for(int i=0; i<num_elements; i++){
                data[i] -= dec_data[i];
                if(fabs(data[i]) > tmp) tmp = fabs(data[i]);
            }
            // std::cout << "max diff = " << tmp << std::endl;
            free(dec_data);
            return approximator_file_size;
        }

        void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            if(filename.size()) approximator_file_name = filename;
        	SZ3::Config conf;
            size_t num = 0;
            auto cmpData = MGARD::readfile<char>(approximator_file_name.c_str(), num);
            approximator_file_size = num;
            SZ_decompress<T>(conf, cmpData.data(), num, data);
        }

        size_t get_size() const {
            return approximator_file_size;
        }

        void print() const {
            std::cout << "PDR dummy approximator" << std::endl;
        }
    private:
        size_t approximator_file_size;
        std::string approximator_file_name = "approximator.dat";
    };
}

#endif