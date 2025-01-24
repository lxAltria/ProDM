#ifndef _PDR_GE_APPROXIMATOR_HPP
#define _PDR_GE_APPROXIMATOR_HPP

#include "ApproximatorInterface.hpp"
#include "SZ3/api/sz.hpp"
#include "utils.hpp"
#include <vector>

namespace PDR {
    // GE approximator with mean-based prediction in vertex-grouped blocks
    template<class T>
    class GEApproximator : public concepts::ApproximatorInterface<T> {
    public:
        GEApproximator(){}
        size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximator_eb, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {

            // block size file
            size_t num_blocks = 0;
            auto block_sizes = MGARD::readfile<int>("/Users/xliang/Github/ProDM_wenbo/GE_small/data/block_sizes.dat", num_blocks);
            std::vector<T> means(num_blocks);
            const T * data_pos = data;
            for(int i=0; i<num_blocks; i++){
                double mean = 0;
                for(int j=0; j<block_sizes[i]; j++){
                    mean += *(data_pos++);
                }
                means[i] = mean / block_sizes[i]; 
            }
            std::cout << "num_blocks = " << num_blocks << std::endl;
            std::cout << "data_pos - data = " << data_pos - data << std::endl;
            std::cout << "dims = " << dimensions[0] << std::endl;
            assert(data_pos - data == dimensions[0]);
            // compress means
            char *cmpData = NULL;
            size_t cmpSize = 0;
            {
                SZ3::Config conf(num_blocks);
                conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
                conf.errorBoundMode = SZ3::EB_REL;
                conf.relErrorBound = 1E-3;
                cmpData = SZ_compress<T>(conf, means.data(), cmpSize);
            }
            // decompress and overwrite
            {
                SZ3::Config conf;
                T * dec_mean = means.data();
                SZ_decompress<T>(conf, cmpData, cmpSize, dec_mean);
            }
            if(filename.size()) approximator_file_name = filename;
            MGARD::writefile(approximator_file_name.c_str(), cmpData, cmpSize);
            approximator_file_size = cmpSize;
            std::cout << "Approximator size = " << approximator_file_size << std::endl;
            T * data_pos2 = data;
            for(int i=0; i<num_blocks; i++){
                for(int j=0; j<block_sizes[i]; j++){
                    *data_pos2 = *data_pos2 - means[i];
                    data_pos2 ++;
                }
            }
            assert(data_pos - data == dimensions[0]);            
            return approximator_file_size;
        }

        void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            if(filename.size()) approximator_file_name = filename;
        	SZ3::Config conf;
            size_t num_blocks = 0;
            auto block_sizes = MGARD::readfile<int>("/Users/xliang/Github/ProDM_wenbo/GE_small/data/block_sizes.dat", num_blocks);
            size_t num = 0;
            auto cmpData = MGARD::readfile<char>(approximator_file_name.c_str(), num);
            approximator_file_size = num;
            std::vector<T> means(num_blocks);
            T * dec_mean = means.data();
            SZ_decompress<T>(conf, cmpData.data(), num, dec_mean);
            T * data_pos = data;
            for(int i=0; i<num_blocks; i++){
                for(int j=0; j<block_sizes[i]; j++){
                    *(data_pos++) = means[i];
                }
            }            
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