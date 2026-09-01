#ifndef _MGARD_SAMPLE_HPP
#define _MGARD_SAMPLE_HPP

#include <vector>
#include <cstdint>

namespace MGARD{

using namespace std;

template<class T>
void profiling_blocks(T const * data, const std::vector<size_t>& dims, std::vector<std::vector<size_t>>& starts, size_t block_size, double releb, size_t stride=4){
    size_t num_dims = dims.size();
    switch (num_dims){
        case 1: {
            size_t dimx = dims[0];
            if(dimx < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            double abseb;
            {
                size_t num_elements = dimx;
                T max = data[0];
                T min = data[0];
                for(int i=1; i<num_elements; i++){
                    if(data[i] > max) max = data[i];
                    if(data[i] < min) min = data[i];
                }
                double abseb = (max - min) * releb;
            }
            for(size_t i = 0; i < dimx - block_size; i += block_size){
                size_t start_idx = i;
                T min = data[start_idx];
                T max = data[start_idx];
                for(size_t ii = 0; ii <= block_size; ii += stride){
                    size_t cur_idx = start_idx + ii;
                    T cur_value = data[cur_idx];
                    if(cur_value < min) min = cur_value;
                    else if(cur_value > max) max = cur_value;
                }
                if(max - min > abseb){
                    size_t a[1] = {i};
                    starts.push_back(std::vector<size_t>(a, a+1));
                }
            }
            break;
        }
        case 2: {
            size_t dimx = dims[0], dimy = dims[1];
            if(dimx < block_size || dimy < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            double abseb;
            {
                size_t num_elements = dimx * dimy;
                T max = data[0];
                T min = data[0];
                for(int i=1; i<num_elements; i++){
                    if(data[i] > max) max = data[i];
                    if(data[i] < min) min = data[i];
                }
                double abseb = (max - min) * releb;
            }
            for(size_t i = 0; i < dimx - block_size; i += block_size){
                for(size_t j = 0; j < dimy - block_size; j += block_size){
                    size_t start_idx = i * dimy + j;
                    T min = data[start_idx];
                    T max = data[start_idx];
                    for(size_t ii = 0; ii <= block_size; ii += stride){
                        for(size_t jj = 0; jj <= block_size; jj += stride){
                            size_t cur_idx = start_idx + ii * dimy + jj;
                            T cur_value = data[cur_idx];
                            if(cur_value < min) min = cur_value;
                            else if(cur_value > max) max = cur_value;
                        }
                    }
                    if(max - min > abseb){
                        size_t a[2] = {i, j};
                        starts.push_back(std::vector<size_t>(a, a+2));
                    }
                }
            }
            break;
        }
        case 3: {
            size_t dimx = dims[0], dimy = dims[1], dimz = dims[2], dimyz = dimy * dimz;
            if(dimx < block_size || dimy < block_size || dimz < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            double abseb;
            {
                size_t num_elements = dimx * dimy * dimz;
                T max = data[0];
                T min = data[0];
                for(int i=1; i<num_elements; i++){
                    if(data[i] > max) max = data[i];
                    if(data[i] < min) min = data[i];
                }
                double abseb = (max - min) * releb;
            }
            for(size_t i = 0; i < dimx - block_size; i += block_size){
                for(size_t j = 0; j < dimy - block_size; j += block_size){
                    for(size_t k = 0; k < dimz - block_size; k += block_size){
                        size_t start_idx = i * dimyz + j * dimz + k;
                        T min = data[start_idx];
                        T max = data[start_idx];
                        for(size_t ii = 0; ii <= block_size; ii += stride){
                            for(size_t jj = 0; jj <= block_size; jj += stride){
                                for(size_t kk = 0; kk <= block_size; kk += stride){
                                    size_t cur_idx = start_idx + ii * dimyz + jj * dimz + kk;
                                    T cur_value = data[cur_idx];
                                    if(cur_value < min) min = cur_value;
                                    else if (cur_value > max) max = cur_value;
                                }
                            }
                        }
                        if(max - min > abseb){
                            size_t a[3] = {i, j, k};
                            starts.push_back(std::vector<size_t>(a, a+3));
                        }
                    }
                }
            }
            break;
        }
        default:
            break;
    }
}

template<class T>
void sample_blocks_after_profiling(T const * data, const std::vector<size_t>& dims, std::vector<std::vector<T>>& sampled_blocks, std::vector<std::vector<size_t>>& starts, size_t block_size, double sample_rate=0.1){
    size_t num_dims = dims.size();
    size_t totalblock_num = 1;
    for(int i = 0; i < num_dims; i++){
        totalblock_num *= static_cast<int>((dims[i] - 1) / block_size);
    }
    size_t num_filtered_blocks = starts.size();
    size_t sample_stride = static_cast<size_t>(num_filtered_blocks / (totalblock_num * sample_rate));
    if (sample_stride <= 0) sample_stride = 1;
    switch(num_dims){
        case 1: {
            for(size_t b = 0; b < num_filtered_blocks; b += sample_stride){
                std::vector<T> sample(block_size, 0);
                size_t startx = starts[b][0];
                for(size_t i = 0; i < block_size; i++){
                    size_t sample_idx = i;
                    size_t idx = i + startx;
                    sample[sample_idx] = data[idx];
                }
                sampled_blocks.push_back(sample);
            }
            break;
        }
        case 2: {
            for(size_t b = 0; b < num_filtered_blocks; b += sample_stride){
                std::vector<T> sample(block_size * block_size, 0);
                size_t startx = starts[b][0], starty = starts[b][1], dimy = dims[1];
                for(size_t i = 0; i < block_size; i++){
                    for(size_t j = 0; j < block_size; j++){
                        size_t sample_idx = i * block_size + j;
                        size_t idx = (i + startx) * dimy + (j + starty);
                        sample[sample_idx] = data[idx];
                    }
                }
                sampled_blocks.push_back(sample);
            }
            break;
        }
        case 3: {
            for(size_t b = 0; b < num_filtered_blocks; b += sample_stride){
                std::vector<T> sample(block_size * block_size * block_size, 0);
                size_t startx = starts[b][0], starty = starts[b][1], startz = starts[b][2], dimy = dims[1], dimz = dims[2];
                size_t square_block_size = block_size * block_size, dimyz = dimy * dimz;
                for(size_t i = 0; i < block_size; i++){
                    for(size_t j = 0; j < block_size; j++){
                        for(size_t k = 0;  k < block_size; k++){
                            size_t sample_idx = i * square_block_size + j * block_size + k;
                            size_t idx = (i + startx) * dimyz + (j + starty) * dimz + k + startz;
                            sample[sample_idx] = data[idx];
                        }
                    }
                }
                sampled_blocks.push_back(sample);
            }
            break;
        }
        default:
            break;
    }
}

template<class T>
void sample_blocks(T const * data, const std::vector<size_t>& dims, std::vector<std::vector<T>>& sampled_blocks, size_t stride=15, size_t block_size=5){
    size_t num_dims = dims.size();
    assert(stride >= 2*block_size);
    switch (num_dims) {
        case 1: {
            if (dims[0] < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            size_t block_1d_size = block_size;
            std::vector<size_t> starts;
            for(size_t i = block_size; i < dims[0] - block_size; i+=stride){
                starts.push_back(i);
            }
            sampled_blocks.resize(starts.size(), std::vector<T>(block_1d_size));
            for(size_t b = 0; b < starts.size(); b++){
                size_t cur_start = starts[b];
                for(size_t i = 0; i < block_size; i++){
                    size_t idx = cur_start + i;
                    sampled_blocks[b][i] = data[idx];
                }
            }
            break;
        }
        case 2: {
            if (dims[0] < block_size || dims[1] < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            size_t stride_n1 = dims[1];
            size_t block_2d_size = block_size * block_size;
            std::vector<size_t> starts;
            for(size_t i = block_size; i < dims[0] - block_size; i+=stride){
                for(size_t j = block_size; j < dims[1] - block_size; j+= stride){
                    starts.push_back(i * stride_n1 + j);
                }
            }
            sampled_blocks.resize(starts.size(), std::vector<T>(block_2d_size));
            for(size_t b = 0 ; b < starts.size(); b++){
                size_t cur_start = starts[b];
                size_t count = 0;
                for(size_t i = 0; i < block_size; i++){
                    for(size_t j = 0; j < block_size; j++){
                        size_t idx = cur_start + i * stride_n1 + j;
                        sampled_blocks[b][count++] = data[idx];
                    }
                }
            }
            break;
        }
        case 3: {
            if (dims[0] < block_size || dims[1] < block_size || dims[2] < block_size){
                std::cerr << "Sampling error: block_size greater than data shape" << std::endl;
                exit(-1);
            }
            size_t stride_n1 = dims[1] * dims[2];
            size_t stride_n2 = dims[2];
            size_t block_3d_size = block_size * block_size * block_size;
            std::vector<size_t> starts;
            for(size_t i = block_size; i < dims[0] - block_size; i+=stride){
                for(size_t j = block_size; j < dims[1] - block_size; j+=stride){
                    for(size_t k = block_size; k < dims[2] - block_size; k+=stride){
                        starts.push_back(i * stride_n1 + j * stride_n2 + k);
                    }
                }
            }
            sampled_blocks.resize(starts.size(), std::vector<T>(block_3d_size));
            for(size_t b = 0; b < starts.size(); b++){
                size_t cur_start = starts[b];
                size_t count = 0;
                // if(b < 2) std::cout << "\nblock " << b << ":\n";
                for(size_t i = 0; i < block_size; i++){
                    for(size_t j = 0; j < block_size; j++){
                        for(size_t k = 0; k < block_size; k++){
                            size_t idx = cur_start + i * stride_n1 + j * stride_n2 + k;
                            sampled_blocks[b][count++] = data[idx];
                            // if (b < 2) std::cout << "index within = \t" << count - 1 << ", real index = " << int(idx) << ", value = " << data[idx] << "\n";
                        }
                    }
                }
                // if(b < 2) std::cout << "\n";
            }
            break;
        }
        default:
            break;
    }
}

}
#endif