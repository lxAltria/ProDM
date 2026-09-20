#ifndef PRODM_DECOMPOSER_MULTILEVEL_MGARDX_UTILS_HPP
#define PRODM_DECOMPOSER_MULTILEVEL_MGARDX_UTILS_HPP

#include <iostream>

#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cassert>
#include <cstdint>

#include "ProDM/Namespace.hpp"
#include "ProDM/Utils/IOUtils.hpp"   // readfile / writefile / print_statistics (used unqualified below)
#include "ProDM/Utils/StatUtils.hpp"

namespace ProDM::MGARDx {

using namespace std;

// compute dimensions for each level
/*
@params dims: dimensions
@params target_level: number of levels to perform
*/
vector<vector<size_t>> init_levels(const vector<size_t>& dims, size_t target_level){
    vector<vector<size_t>> level_dims;
    // compute n_nodal in each level
    for(int i=0; i<=target_level; i++){
        level_dims.push_back(vector<size_t>(dims.size()));
    }
    for(int i=0; i<dims.size(); i++){
        int n = dims[i];
        for(int j=0; j<=target_level; j++){
            level_dims[target_level - j][i] = n;
            n = (n >> 1) + 1;
        }
    }
    return level_dims;
}
template <class T>
inline T interp_cubic(T a, T b, T c, T d) {
    return (-a + 9 * b + 9 * c - d) / 16;
}
template <class T>
inline T interp_quad_1(T a, T b, T c) {
    return (3 * a + 6 * b - c) / 8;
}
template <class T>
inline T interp_quad_2(T a, T b, T c) {
    return (-a + 6 * b + 3 * c) / 8;
}
template<class T>
inline T interp_quad_3(T a, T b, T c) {
    // return (3 * a - 10 * b + 15 * c) / 8;
    return (3 * a + 6 * b - c) / 8;
}
std::vector<std::vector<uint32_t>> compute_level_dims_new(const std::vector<uint32_t>& dims, uint32_t target_level){
    std::vector<std::vector<uint32_t>> level_dims;
    for(int i=0; i<=target_level; i++){
        level_dims.push_back(std::vector<uint32_t>(dims.size()));
    }
    for(int i=0; i<dims.size(); i++){
        int n = dims[i];
        for(int j=0; j<=target_level; j++){
            level_dims[target_level - j][i] = n;
            n = (n >> 1) + (n & 1);
        }
    }
    return level_dims;
}
std::vector<uint32_t> compute_level_buffers_size_2D_coeff(const std::vector<std::vector<uint32_t>>& level_dims, int target_level, std::vector<std::vector<uint32_t>>& level_buffer_dims){
    assert(level_dims.size());
    size_t num_dims = level_dims[0].size();
    // size_t count;
    size_t size;
    std::vector<uint32_t> level_sizes(1 + num_dims * target_level);
    level_sizes[0] = 1;
    std::vector<uint32_t> temp_level_dims;
    for(int i=0; i<num_dims; i++){
        level_sizes[0] *= level_dims[0][i];
        temp_level_dims.push_back(level_dims[0][i]);
    }
    level_buffer_dims.push_back(temp_level_dims);
    for(size_t l=1; l<=target_level; l++){
        // count = 0;
        switch (num_dims){
            case 1:
            {
                size = level_dims[l][0];
                temp_level_dims = level_dims[l];
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                // assert(count == (level_dims[l][0] - level_dims[l-1][0]));
                break;
            }
            case 2:
            {
                size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1];
                temp_level_dims = {(level_dims[l][0] - level_dims[l-1][0]), level_dims[l-1][1]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]);
                temp_level_dims = {level_dims[l][0], (level_dims[l][1] - level_dims[l-1][1])};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 2] = size;
                // count += size;
                // assert(count == (level_dims[l][0]*level_dims[l][1] - level_dims[l-1][0]*level_dims[l-1][1]));
                break;
            }
            case 3:
            {
                // n1 (cur_n1 - pre_n1) * pre_n2 * pre_n3
                size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1] * level_dims[l-1][2];
                temp_level_dims = {(level_dims[l][0] - level_dims[l-1][0]), level_dims[l-1][1], level_dims[l-1][2]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                // n2 cur_n1 * (cur_n2 - pre_n2) * pre_n3
                size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]) * level_dims[l-1][2];
                temp_level_dims = {level_dims[l][0], (level_dims[l][1] - level_dims[l-1][1]), level_dims[l-1][2]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 2] = size;
                // count += size;
                // n3 cur_n1 * cur_n2 * (cur_n3 - pre_n3)
                size = level_dims[l][0] * level_dims[l][1] * (level_dims[l][2] - level_dims[l-1][2]);
                temp_level_dims = {level_dims[l][0], level_dims[l][1], (level_dims[l][2] - level_dims[l-1][2])};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 3] = size;
                // count += size;
                // assert(count == (level_dims[l][0]*level_dims[l][1]*level_dims[l][2] - level_dims[l-1][0]*level_dims[l-1][1]*level_dims[l-1][2]));
                break;
            }
            default:
                std::cerr << num_dims << "-Dimentional decomposition not implemented." << std::endl;
                exit(-1);
        }
    }
    return level_sizes;
}
std::vector<uint32_t> compute_level_buffers_size(const std::vector<std::vector<uint32_t>>& level_dims, int target_level, std::vector<std::vector<uint32_t>>& level_buffer_dims){
    assert(level_dims.size());
    size_t num_dims = level_dims[0].size();
    // size_t count;
    size_t size;
    std::vector<uint32_t> level_sizes(1 + num_dims * target_level);
    level_sizes[0] = 1;
    std::vector<uint32_t> temp_level_dims;
    for(int i=0; i<num_dims; i++){
        level_sizes[0] *= level_dims[0][i];
        temp_level_dims.push_back(level_dims[0][i]);
    }
    level_buffer_dims.push_back(temp_level_dims);
    for(size_t l=1; l<=target_level; l++){
        // count = 0;
        switch (num_dims){
            case 1:
            {
                size = level_dims[l][0];
                temp_level_dims = level_dims[l];
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                // assert(count == (level_dims[l][0] - level_dims[l-1][0]));
                break;
            }
            case 2:
            {
                size = level_dims[l-1][0] * (level_dims[l][1] - level_dims[l-1][1]);
                temp_level_dims = {level_dims[l-1][0], (level_dims[l][1] - level_dims[l-1][1])};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l][1];
                temp_level_dims = {(level_dims[l][0] - level_dims[l-1][0]), level_dims[l][1]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 2] = size;
                // count += size;
                // assert(count == (level_dims[l][0]*level_dims[l][1] - level_dims[l-1][0]*level_dims[l-1][1]));
                break;
            }
            case 3:
            {
                // interp direction: n3, dims: pre_n1 * pre_n2 * (cur_n3 - pre_n3)
                size = level_dims[l-1][0] * level_dims[l-1][1] * (level_dims[l][2] - level_dims[l-1][2]);
                temp_level_dims = {level_dims[l-1][0], level_dims[l-1][1], (level_dims[l][2] - level_dims[l-1][2])};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 1] = size;
                // count += size;
                // interp direction: n2, dims: pre_n1 * (cur_n2 - pre_n2) * cur_n3
                size = level_dims[l-1][0] * (level_dims[l][1] - level_dims[l-1][1]) * level_dims[l][2];
                temp_level_dims = {level_dims[l-1][0], (level_dims[l][1] - level_dims[l-1][1]), level_dims[l][2]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 2] = size;
                // count += size;
                // interp direction: n1, dims: (cur_n3 - pre_n3) * cur_n2 * cur_n3
                size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l][1] * level_dims[l][2];
                temp_level_dims = {(level_dims[l][0] - level_dims[l-1][0]), level_dims[l][1], level_dims[l][2]};
                level_buffer_dims.push_back(temp_level_dims);
                level_sizes[(l-1) * num_dims + 3] = size;
                // count += size;
                // assert(count == (level_dims[l][0]*level_dims[l][1]*level_dims[l][2] - level_dims[l-1][0]*level_dims[l-1][1]*level_dims[l-1][2]));
                break;
            }
            default:
                std::cerr << num_dims << "-Dimentional decomposition not implemented." << std::endl;
                exit(-1);
        }
    }
    return level_sizes;
}
// Generic buffer size computation for arbitrary interpolation direction order.
// interp_order: permutation of {0, 1, ..., num_dims-1} specifying the interpolation order.
//   e.g. {0,1,2} = X->Y->Z,  {2,1,0} = Z->Y->X
// For sub-level s (interpolating dimension interp_order[s]):
//   - dimensions already interpolated (interp_order[0..s-1]): use cur resolution
//   - dimension being interpolated (interp_order[s]): use cur - pre (new points only)
//   - dimensions not yet interpolated (interp_order[s+1..]): use pre resolution
std::vector<uint32_t> compute_level_buffers_size_generic(
    const std::vector<std::vector<uint32_t>>& level_dims, int target_level,
    const std::vector<uint32_t>& interp_order,
    std::vector<std::vector<uint32_t>>& level_buffer_dims){
    assert(level_dims.size());
    size_t num_dims = level_dims[0].size();
    assert(interp_order.size() == num_dims);
    std::vector<uint32_t> level_sizes(1 + num_dims * target_level);
    level_sizes[0] = 1;
    std::vector<uint32_t> temp_level_dims(num_dims);
    for(size_t i=0; i<num_dims; i++){
        level_sizes[0] *= level_dims[0][i];
        temp_level_dims[i] = level_dims[0][i];
    }
    level_buffer_dims.push_back(temp_level_dims);
    for(size_t l=1; l<=target_level; l++){
        std::vector<bool> interpolated(num_dims, false);
        for(size_t s=0; s<num_dims; s++){
            uint32_t interp_dim = interp_order[s];
            size_t size = 1;
            for(size_t d=0; d<num_dims; d++){
                if(d == interp_dim){
                    temp_level_dims[d] = level_dims[l][d] - level_dims[l-1][d];
                } else if(interpolated[d]){
                    temp_level_dims[d] = level_dims[l][d];
                } else {
                    temp_level_dims[d] = level_dims[l-1][d];
                }
                size *= temp_level_dims[d];
            }
            level_buffer_dims.push_back(temp_level_dims);
            level_sizes[(l-1) * num_dims + s + 1] = size;
            interpolated[interp_dim] = true;
        }
    }
    return level_sizes;
}
}
#endif