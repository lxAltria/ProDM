#ifndef _PDR_MGARD_APPROXIMATOR_HPP
#define _PDR_MGARD_APPROXIMATOR_HPP

#undef REL
#undef ABS

#include "ApproximatorInterface.hpp"
#include "mgard/compress_x.hpp"
#include "utils.hpp"

namespace PDR {
    // MGARD approximator with MGARDx prediction
    template<class T>
    class MGARDApproximator : public concepts::ApproximatorInterface<T> {
    public:
        MGARDApproximator(){}
        size_t refactor_approximate(T * data, const std::vector<uint32_t>& dimensions, T approximator_eb, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            size_t cmpSize = 0;
            void *cmpData = nullptr;
            size_t num_elements = std::accumulate(dimensions.begin(), dimensions.end(), size_t{1}, std::multiplies<size_t>());
            mgard_x::SIZE n1 = 0;
            mgard_x::SIZE n2 = 0;
            mgard_x::SIZE n3 = 0;
            mgard_x::data_type dtype = std::is_same_v<T, double> ? mgard_x::data_type::Double : mgard_x::data_type::Float;
            mgard_x::Config config;
            config.lossless = mgard_x::lossless_type::Huffman_Zstd;
            config.dev_type = mgard_x::device_type::SERIAL;
            double s = 0, norm;
            if(dimensions.size() == 1){
                n1 = dimensions[0];
                std::vector<mgard_x::SIZE> shape{n1};
                mgard_x::compress(1, dtype, shape, approximator_eb, s, mgard_x::error_bound_type::ABS, data, cmpData, cmpSize, config, false);
            }
            else if(dimensions.size() == 2){
                n1 = dimensions[0];
                n2 = dimensions[1];
                std::vector<mgard_x::SIZE> shape{n1, n2};
                mgard_x::compress(2, dtype, shape, approximator_eb, s, mgard_x::error_bound_type::ABS, data, cmpData, cmpSize, config, false);
            }
            else if(dimensions.size() == 3){
                n1 = dimensions[0];
                n2 = dimensions[1];
                n3 = dimensions[2];
                std::vector<mgard_x::SIZE> shape{n1, n2, n3};
                mgard_x::compress(3, dtype, shape, approximator_eb, s, mgard_x::error_bound_type::ABS, data, cmpData, cmpSize, config, false);
            }
            else{
                std::cout << "Dimension larger than 4 is not supported at SZ approximator!" << std::endl;
                exit(-1);
            }
            if(filename.size()) approximator_file_name = filename;
            MGARD::writefile(approximator_file_name.c_str(), static_cast<char*>(cmpData), cmpSize);
            approximator_file_size = cmpSize;
            // std::cout << "Approximator size = " << approximator_file_size << std::endl;
            // std::cout << "num_elements = " << num_elements << std::endl;
            void * void_dec_data = nullptr;
            mgard_x::decompress(cmpData, cmpSize, void_dec_data, false);
            T * dec_data = static_cast<T*>(void_dec_data);
            free(cmpData);
            // T tmp = 0;
            for(int i=0; i<num_elements; i++){
                data[i] -= dec_data[i];
            //     if(fabs(data[i]) > tmp) tmp = fabs(data[i]);
            }
            // std::cout << "max diff = " << tmp << std::endl;
            free(dec_data);
            return approximator_file_size;
        }

        void reconstruct_approximate(T * data, const std::vector<uint32_t>& dimensions, std::string filename=std::string(""), std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            if(filename.size()) approximator_file_name = filename;
        	mgard_x::Config config;
            config.lossless = mgard_x::lossless_type::Huffman_Zstd;
            config.dev_type = mgard_x::device_type::SERIAL;
            size_t num = 0;
            auto cmpData = MGARD::readfile<char>(approximator_file_name.c_str(), num);
            approximator_file_size = num;

            void * void_dec_data = nullptr;
            mgard_x::decompress(cmpData.data(), num, void_dec_data, config, false);

            size_t num_elements = std::accumulate(dimensions.begin(), dimensions.end(), size_t{1}, std::multiplies<size_t>());
            memcpy(data, void_dec_data, num_elements * sizeof(T));
            free(void_dec_data);
        }

        size_t get_size() const {
            return approximator_file_size;
        }

        void print() const {
            std::cout << "PDR MGARD approximator" << std::endl;
        }
    private:
        size_t approximator_file_size;
        std::string approximator_file_name = "approximator.dat";
    };
}

#endif