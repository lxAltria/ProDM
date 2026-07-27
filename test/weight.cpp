#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "/Users/wenboli/uky/ProDM/external/MGARDx/include/utils.hpp"
#include "/Users/wenboli/uky/ProDM/external/MGARDx/include/reorder.hpp"
#include "MDR/RefactorUtils.hpp"
#include <unistd.h>

using namespace std;
template <class T>
void switch_rows_2D_by_buffer(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    T * nodal_data_buffer = data_buffer + n2; // skip the first nodal row
    T * coeff_data_buffer = data_buffer + n1_nodal * n2;
    T * cur_data_pos = data_pos + stride;
    for(int i=0; i<n1_coeff; i++){
        // copy coefficient rows
        memcpy(coeff_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
        cur_data_pos += stride;
        // copy nodal rows
        memcpy(nodal_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    if(!(n1&1)){
        // n1 is even, move the last nodal row
        memcpy(coeff_data_buffer - n2, cur_data_pos, n2 * sizeof(T));
    }
    // copy data back
    cur_data_pos = data_pos + stride;
    for(int i=1; i<n1; i++){
        memcpy(cur_data_pos, data_buffer + i * n2, n2 * sizeof(T));
        cur_data_pos += stride;
    }
}
template <class T>
void weight_reorder_1D(const T * data_pos, size_t n_nodal, size_t n_coeff, T * nodal_buffer, T * coeff_buffer){
    T * nodal_pos = nodal_buffer;
    T * coeff_pos = coeff_buffer;
    T const * cur_data_pos = data_pos;
    for(int i=0; i<n_coeff; i++){
        *(nodal_pos++) = *(cur_data_pos++);
        *(coeff_pos++) = *(cur_data_pos++);
    }
    *(nodal_pos++) = *(cur_data_pos++);
    if(n_nodal == n_coeff + 2){
        // *nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
        // modified process
        *nodal_pos = *cur_data_pos;
    }
}
/*
    oxoxo       oooxx       oooxx
    xxxxx   (1) xxxxx   (2) oooxx
    oxoxo   =>  oooxx   =>  oooxx
    xxxxx       xxxxx       xxxxx
    oxoxo       oooxx       xxxxx
*/
template <class T>
void weight_reorder_2D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    for(int i=0; i<n1; i++){
        weight_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    // do reorder (2)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
}
template <class T>
void weight_reorder_3D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // do 2D reorder
    for(int i=0; i<n1; i++){
        weight_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
}
template <class T>
inline void compare_and_restrict(T& a, T b){
    if(a > b) a = b;
}
template <class T>
void compare_and_restrict_1D(int n, T * weight){
    // validate weight[0]
    std::cout << n << std::endl;
    compare_and_restrict(weight[0], weight[1]);
    // validate all nodal nodes
    for(int i=2; i<n-1; i+=2){
        // if(i == 4) std::cout << weight[i-1] << " " << weight[i] << " " << weight[i+1] << std::endl;
        compare_and_restrict(weight[i], weight[i-1]);
        compare_and_restrict(weight[i], weight[i+1]);
    }
    if(n & 1){
        // validate the last node
        compare_and_restrict(weight[n-1], weight[n-2]);
    }
}
template <class T>
void compare_and_restrict_2D(int n1, int n2, size_t stride, T * weight){
    T * weight_pos = weight;
    // validate first row
    compare_and_restrict(weight_pos[0], weight_pos[1]);
    compare_and_restrict(weight_pos[0], weight_pos[stride]);
    compare_and_restrict(weight_pos[0], weight_pos[stride+1]);
    for(int j=2; j<n2-1; j+=2){
        compare_and_restrict(weight_pos[j], weight_pos[j-1]);
        compare_and_restrict(weight_pos[j], weight_pos[j+1]);           
        compare_and_restrict(weight_pos[j], weight_pos[j-1+stride]);
        compare_and_restrict(weight_pos[j], weight_pos[j+1+stride]);
    }
    if(n2 & 1){
        compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2]);
        compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2 + stride]);
        compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 1 + stride]);       
    }
    weight_pos += 2*stride;
    // validate all nodal rows
    for(int i=2; i<n1-1; i+=2){
        compare_and_restrict(weight_pos[0], weight_pos[-stride]);
        compare_and_restrict(weight_pos[0], weight_pos[-stride+1]);
        compare_and_restrict(weight_pos[0], weight_pos[stride]);
        compare_and_restrict(weight_pos[0], weight_pos[stride+1]);
        for(int j=2; j<n2-1; j+=2){
            compare_and_restrict(weight_pos[j], weight_pos[j-stride-1]);
            compare_and_restrict(weight_pos[j], weight_pos[j-stride+1]);
            compare_and_restrict(weight_pos[j], weight_pos[j+stride-1]);
            compare_and_restrict(weight_pos[j], weight_pos[j+stride+1]);
        }
        if(n2 & 1){
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2 - stride]);
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 1 - stride]);
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2 + stride]);
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 1 + stride]);
        }
        weight_pos += 2*stride;
    }
    if(n1 & 1){
        compare_and_restrict(weight_pos[0], weight_pos[-stride]);
        compare_and_restrict(weight_pos[0], weight_pos[-stride+1]);
        compare_and_restrict(weight_pos[0], weight_pos[1]);
        for(int j=2; j<n2-1; j+=2){
            compare_and_restrict(weight_pos[j], weight_pos[j-stride-1]);
            compare_and_restrict(weight_pos[j], weight_pos[j-stride+1]);
            compare_and_restrict(weight_pos[j], weight_pos[j-1]);
            compare_and_restrict(weight_pos[j], weight_pos[j+1]);           
        }
        if(n2 & 1){
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 1 - stride]);
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2 - stride]);
            compare_and_restrict(weight_pos[n2 - 1], weight_pos[n2 - 2]);
        }       
    }
}
template <class T>
void propagateWeightInLevel(const vector<uint32_t>& dims, const vector<uint32_t>& strides, vector<T>& weight){
    int num_dims = dims.size();
    if(num_dims == 1){
        compare_and_restrict_1D(dims[0], weight.data());
    }
    else if(num_dims == 2){
        compare_and_restrict_2D(dims[0], dims[1], strides[1], weight.data());
    }
    else if(num_dims == 3){
        //compare_and_restrict_3D(dims[0], dims[1], dims[2], strides[0], strides[1], weight.data());
    }
    else{
        std::cout << "dimension high than 4 is not supported\n";
    }
}
template <class T>
void propagateWeight(const vector<uint32_t>& dims, int target_level, vector<T>& weight){
    int num_dims = dims.size();
    vector<uint32_t> strides(num_dims);
    uint32_t stride = 1;
    for(int i=0; i<num_dims; i++){
        strides[num_dims-1-i] = stride;
        stride *= dims[num_dims-1-i];
    }
    vector<T> buffer(stride);
    T * data_buffer = buffer.data();
    T * data_pos = weight.data();
    auto level_dims = MDR::compute_level_dims(dims, target_level);
    // MDR::print_vec("level_dims", level_dims);
    std::cout << "start propagation\n";
    for(int i=0; i<target_level; i++){
        // MDR::print_vec(weight);
        propagateWeightInLevel(level_dims[target_level - i], strides, weight);
        // MDR::print_vec(weight);
        if(num_dims == 1){
            size_t n = level_dims[target_level - i][0];
            size_t n_nodal = (n >> 1) + 1;
            size_t n_coeff = n - n_nodal;
            T * nodal_buffer = data_buffer;
            T * coeff_buffer = data_buffer + n_nodal;
            weight_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
            memcpy(data_pos, data_buffer, n*sizeof(T));           
        }
        else if (num_dims == 2){
            size_t n1 = level_dims[target_level - i][0];
            size_t n2 = level_dims[target_level - i][1];
            weight_reorder_2D(data_pos, data_buffer, n1, n2, strides[0]);
        }
        else if(num_dims == 3){
            size_t n1 = level_dims[target_level - i][0];
            size_t n2 = level_dims[target_level - i][1];
            size_t n3 = level_dims[target_level - i][2];
            // weight_reorder_3D(data_pos, data_buffer, n1, n2, n3, strides[0], strides[1]);
        }
        else{
            std::cout << "dimension high than 4 is not supported\n";
        }
        // MDR::print_vec(weight);
        std::cout << "level done\n";
    }
}
template <class T>
void assign_block_value_1D(const size_t n1, const int block_size, T * data){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    T * data_x_pos = data;
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        T * cur_data_pos = data_x_pos;
        // find block maximum
        T block_max = *cur_data_pos;
        for(int ii=0; ii<size_1; ii++){
            if(*cur_data_pos > block_max) block_max = *cur_data_pos;
            cur_data_pos ++;
        }
        // assign block maximum
        cur_data_pos = data_x_pos;
        for(int ii=0; ii<size_1; ii++){
            *(cur_data_pos ++) = block_max;
        }
        data_x_pos += size_1;
    }
}
template <class T>
void assign_block_value_2D(const size_t n1, const size_t n2, const size_t dim0_offset, const int block_size, T * data){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    T * data_x_pos = data;
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        T * data_y_pos = data_x_pos;
        for(int j=0; j<num_block_2; j++){
            int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
            T * cur_data_pos = data_y_pos;
            // find block maximum
            T block_max = *cur_data_pos;
            for(int ii=0; ii<size_1; ii++){
                for(int jj=0; jj<size_2; jj++){
                    if(*cur_data_pos > block_max) block_max = *cur_data_pos;
                    cur_data_pos ++;
                }
                cur_data_pos += dim0_offset - size_2;
            }
            // assign block maximum
            cur_data_pos = data_y_pos;
            for(int ii=0; ii<size_1; ii++){
                for(int jj=0; jj<size_2; jj++){
                    *(cur_data_pos ++) = block_max;
                }
                cur_data_pos += dim0_offset - size_2;
            }
            data_y_pos += size_2;
        }
        data_x_pos += dim0_offset * size_1;
    }
}
template <class T>
void assign_block_value_3D(const size_t n1, const size_t n2, const size_t n3, const size_t dim0_offset, const size_t dim1_offset, const int block_size, T * data){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    size_t num_block_3 = (n3 - 1) / block_size + 1;
    T * data_x_pos = data;
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        T * data_y_pos = data_x_pos;
        for(int j=0; j<num_block_2; j++){
            int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
            T * data_z_pos = data_y_pos;
            for(int k=0; k<num_block_3; k++){
                int size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                T * cur_data_pos = data_z_pos;
                // find block maximum
                T block_max = *cur_data_pos;
                for(int ii=0; ii<size_1; ii++){
                    for(int jj=0; jj<size_2; jj++){
                        for(int kk=0; kk<size_3; kk++){
                            if(*cur_data_pos > block_max) block_max = *cur_data_pos;
                            cur_data_pos ++;
                        }
                        cur_data_pos += dim1_offset - size_3;
                    }
                    cur_data_pos += dim0_offset - size_2 * dim1_offset;
                }
                cur_data_pos = data_z_pos;
                for(int ii=0; ii<size_1; ii++){
                    for(int jj=0; jj<size_2; jj++){
                        for(int kk=0; kk<size_3; kk++){
                            *(cur_data_pos ++) = block_max;
                        }
                        cur_data_pos += dim1_offset - size_3;
                    }
                    cur_data_pos += dim0_offset - size_2 * dim1_offset;
                }
                // assign block maximum
                data_z_pos += size_3;
            }
            data_y_pos += dim1_offset * size_2;
        }
        data_x_pos += dim0_offset * size_1;
    }
}
template <class T>
std::vector<int> normalize_weights(std::vector<T>& weights, int num_weight_bitplane=0){
    int max_weight = 4;
    std::vector<int> int_weights(weights.size());
    T max = fabs(weights[0]);
    T min;
    if (weights[0] == 0){
        for(int i=1; i<weights.size(); i++){
            if(weights[i]){
                min = fabs(weights[i]);
            }
        }
    }
    else{
        min = fabs(weights[0]);
    }
    // std::cout << "weights[0] = " << weights[0] << std::endl; 
    for(int i=1; i<weights.size(); i++){
        if(weights[i]){
            auto abs_weight = fabs(weights[i]);
            if(abs_weight > max) max = abs_weight;
            if(abs_weight < min) min = abs_weight;
        }
    }
    std::cout << max << " " << min << std::endl;
    // exit(0);
    for(int i=0; i<weights.size(); i++){
        if(!weights[i]){
            int_weights[i] = 0; ///// 0
        }
        else {
            bool sign = weights[i] < 0;
            // weights[i] = sign ? -log2(-weights[i]/min) : log2(weights[i]/min);
            int_weights[i] = sign ? static_cast<int>(log2(-weights[i]/min)) : static_cast<int>(log2(weights[i]/min));
            if (int_weights[i] > max_weight) int_weights[i] = max_weight;
        }
        // int_weights[i] = max_weight - int_weights[i];    // for value and QoI, comment it when dQoI
    }
    return int_weights;
}
std::vector<int> get_block_weight_1D(const size_t n1, std::vector<int>& weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    std::vector<int> int_weights(num_block_1);
    for(int i=0; i<num_block_1; i++){
        int_weights[i] = weights[i * block_size];
    }
    return int_weights;
}
std::vector<int> get_block_weight_2D(const size_t n1, const size_t n2, std::vector<int>& weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    std::vector<int> int_weights(num_block_1 * num_block_2);
    for(int i=0; i<num_block_1; i++){
        for (int j=0; j<num_block_2; j++){
            int_weights[i * num_block_2 + j] = weights[i * n2 * block_size + j * block_size];
        }
    }
    return int_weights;
}
std::vector<int> get_block_weight_3D(const size_t n1, const size_t n2, const size_t n3, std::vector<int>& weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    size_t num_block_3 = (n3 - 1) / block_size + 1;
    std::vector<int> int_weights(num_block_1 * num_block_2 * num_block_3);
    for(int i=0; i<num_block_1; i++){
        for(int j=0; j<num_block_2; j++){
            for(int k=0; k<num_block_3; k++){
                int_weights[i * num_block_2 * num_block_3 + j * num_block_3 + k] = weights[i * n2 * n3 * block_size + j * n3 * block_size + k * block_size];
            }
        }
    }
    return int_weights;
}
std::vector<int> fill_block_weight_1D(const size_t n1, std::vector<int> weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    std::vector<int> filled_weight(n1);
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        for(int ii=0; ii<size_1; ii++){
            filled_weight[i * block_size + ii] = weights[i];
        }
    }
    return filled_weight;
}
std::vector<int> fill_block_weight_2D(const size_t n1, const size_t n2, std::vector<int>& weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    std::vector<int> filled_weight(n1 * n2);
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        for(int j=0; j<num_block_2; j++){
            int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
            for(int ii=0; ii < size_1; ii++){
                for(int jj=0; jj < size_2; jj++){
                    filled_weight[i * n2 * block_size + j * block_size + ii * n2 + jj] = weights[i * num_block_2 + j];
                }
            }
        }
    }
    return filled_weight;
}
std::vector<int> fill_block_weight_3D(const size_t n1, const size_t n2, const size_t n3, std::vector<int>& weights, int block_size){
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    size_t num_block_3 = (n3 - 1) / block_size + 1;
    std::vector<int> filled_weight(n1 * n2 * n3);
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        for(int j=0; j<num_block_2; j++){
            int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
            for(int k=0; k<num_block_3; k++){
                int size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                for(int ii=0; ii < size_1; ii++){
                    for(int jj=0; jj < size_2; jj++){
                        for(int kk=0; kk < size_3; kk++){
                            filled_weight[i * n2 * n3 * block_size + j * n3 * block_size + k * block_size + ii * n2 * n3 + jj * n3 + kk] = weights[i * num_block_2 * num_block_3 + j * num_block_3 + k];
                        }
                    }
                }
            }
        }
    }
    return filled_weight;
}
bool compare_weight(const size_t n, std::vector<int> ori_weight, std::vector<int> filled_weight){
    bool same = true;
    int count = 0;
    for(int i=0; i<n; i++){
        if (ori_weight[i] != filled_weight[i]){
            if(same) same = false;
            //std::cout << "ori_weight[" << i << "] = " << ori_weight[i] << " filled_weight[" << i << "] = " << filled_weight[i] << std::endl;
            count++;
        }
    }
    std::cout << "different count = " << count << std::endl;
    return same;
}
void write_weight(const std::string& file_path, const std::vector<int>& block_weights, int block_size) {
    // 输出最大和最小权重以便调试
    int max_w = block_weights[0];
    int min_w = block_weights[0];
    for (int w : block_weights) {
        if (w > max_w) max_w = w;
        if (w < min_w) min_w = w;
    }
    std::cout << "Min Weight: " << min_w << ", Max Weight: " << max_w << std::endl;

    // 计算权重的总大小
    uint32_t weight_size = sizeof(int) + sizeof(size_t) + block_weights.size() * sizeof(int);

    // 分配缓冲区
    uint8_t* weight_data = (uint8_t*)malloc(weight_size);
    uint8_t* weight_data_pos = weight_data;

    // 写入 block_size
    memcpy(weight_data_pos, &block_size, sizeof(int));
    weight_data_pos += sizeof(int);

    // 写入权重数组的大小
    size_t block_weight_size = block_weights.size();
    memcpy(weight_data_pos, &block_weight_size, sizeof(size_t));
    weight_data_pos += sizeof(size_t);

    // 写入权重数组数据
    memcpy(weight_data_pos, block_weights.data(), block_weight_size * sizeof(int));

    // 写入到文件
    FILE* file = fopen(file_path.c_str(), "wb");
    if (!file) {
        perror("Error opening file for writing");
        free(weight_data);
        return;
    }
    fwrite(weight_data, 1, weight_size, file);
    fclose(file);

    free(weight_data);
    std::cout << "Weight data written to " << file_path << std::endl;
}

void load_weight(const std::string& file_path, std::vector<int>& block_weights, int& block_size) {
    // 打开文件读取
    FILE* file = fopen(file_path.c_str(), "rb");
    if (!file) {
        perror("Error opening file for reading");
        return;
    }

    // 获取文件大小
    fseek(file, 0, SEEK_END);
    uint32_t num_bytes = ftell(file);
    rewind(file);

    // 分配缓冲区
    uint8_t* weight_data = (uint8_t*)malloc(num_bytes);
    fread(weight_data, 1, num_bytes, file);
    fclose(file);

    uint8_t* weight_data_pos = weight_data;

    // 读取 block_size
    memcpy(&block_size, weight_data_pos, sizeof(int));
    weight_data_pos += sizeof(int);

    // 读取权重数组的大小
    size_t block_weight_size;
    memcpy(&block_weight_size, weight_data_pos, sizeof(size_t));
    weight_data_pos += sizeof(size_t);

    // 读取权重数组
    block_weights.clear();
    block_weights.assign(reinterpret_cast<int*>(weight_data_pos), reinterpret_cast<int*>(weight_data_pos) + block_weight_size);

    free(weight_data);

    // 输出最大和最小权重以便调试
    int max_w = block_weights[0];
    int min_w = block_weights[0];
    for (int w : block_weights) {
        if (w > max_w) max_w = w;
        if (w < min_w) min_w = w;
    }

}
int main(int argc, char ** argv){
    using T = float;
    using T_stream = uint32_t;
    int arg = 1;
    std::string filename(argv[arg++]);
    int num_dims = atoi(argv[arg++]);
    std::vector<uint32_t> dims(num_dims);
    int n = 1;
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[arg++]);
        n *= dims[i];
    }
    int block_size = atoi(argv[arg++]);
    size_t num_elements = 0;
    auto data = MGARD::readfile<float>(filename.c_str(), num_elements);
    assert(n == num_elements);
    // MDR::print_vec(data);
    if(dims.size() == 1){
        assign_block_value_1D(dims[0], block_size, data.data());
    }
    else if(dims.size() == 2){
        assign_block_value_2D(dims[0], dims[1], dims[1], block_size, data.data());
    }
    else if(dims.size() == 3){
        assign_block_value_3D(dims[0], dims[1], dims[2], dims[1]*dims[2], dims[2], block_size, data.data());
    }
    // MDR::print_vec(data);
    std::vector<int> weight = normalize_weights(data);
    std::vector<int> int_weight;
    if(dims.size() == 1){
        int_weight = get_block_weight_1D(dims[0], weight, block_size);
    }
    else if(dims.size() == 2){
        int_weight = get_block_weight_2D(dims[0], dims[1], weight, block_size);
    }
    else if(dims.size() == 3){
        int_weight = get_block_weight_3D(dims[0], dims[1], dims[2], weight, block_size);
    }
    std::cout << "int_weight.size() = " << int_weight.size() << std::endl;
    std::string file_path = "weight.bin";
    write_weight(file_path, int_weight, block_size);
    std::vector<int> loaded_weights;
    int loaded_block_size;
    load_weight(file_path, loaded_weights, loaded_block_size);
    if(compare_weight(loaded_weights.size(), int_weight, loaded_weights)){
        std::cout << "Loaded Weight Same.\n";
    }
    else{
        std::cout << "Loaded Weight Not the same.\n";
    }
    std::vector<int> filled_weight(n);
    if(dims.size() == 1){
        // filled_weight = fill_block_weight_1D(dims[0], int_weight, block_size);
        filled_weight = fill_block_weight_1D(dims[0], loaded_weights, loaded_block_size);
    }
    else if(dims.size() == 2){
        // filled_weight = fill_block_weight_2D(dims[0], dims[1], int_weight, block_size);
        filled_weight = fill_block_weight_2D(dims[0], dims[1], loaded_weights, loaded_block_size);
    }
    else if(dims.size() == 3){
        // filled_weight = fill_block_weight_3D(dims[0], dims[1], dims[2], int_weight, block_size);
        filled_weight = fill_block_weight_3D(dims[0], dims[1], dims[2], loaded_weights, loaded_block_size);
    }
    // MDR::print_vec(data);
    // propagateWeight(dims, 3, weight);
    // propagateWeight(dims, 3, filled_weight);
    if(compare_weight(n, weight, filled_weight)){
        std::cout << "Filled Weight Same.\n";
    }
    else{
        std::cout << "Filled Weight Not the same.\n";
    }
    // MDR::print_vec(data);
    return 0;
}