#ifndef _MGARD_COEFF_DECOMPOSE_INTERLEAVE_HPP
#define _MGARD_COEFF_DECOMPOSE_INTERLEAVE_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "reorder.hpp"
#include "ProDM/MGARDx/utils.hpp"
#include "correction.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Coeff_Decomposer_Interleaver{
public:
	Coeff_Decomposer_Interleaver(bool use_sz_=true){
            use_sz = use_sz_;
        };
	~Coeff_Decomposer_Interleaver(){
	};
    // return decomposed and interleaved data buffers for each level, num_levels = Dims * MGARD_level (target_level)
	// Combining MGARD Decomposer and Interleaver
	std::vector<std::vector<T>> decompose(T * data, const vector<size_t>& dims, const size_t direction=0, size_t target_level=1, bool hierarchical=false, bool cubic=false, vector<size_t> strides=vector<size_t>()){
		if(dims.size() != 3) {
			std::cerr << "Only support 3D dataset" << std::endl;
			exit(-1);
		}
		data_begin = data;
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
        int max_level = log2(*min_element(dims.begin(), dims.end()));
        if(target_level > max_level) target_level = max_level;
		init(dims, direction, target_level);
		
		size_t h = 1;
		size_t n1 = dims[0];
		size_t n2 = dims[1];
		size_t n3 = dims[2];
		for(int current_level=target_level; current_level >= 0; current_level--){
			decompose_interleave_3D_HB_with_direction(data, n1, n2, n3, h, direction, current_level);
			h <<= 1;
		}
		// std::cout << "avg abs error = " << sum / (num_elements - level_sizes[0]) << std::endl;
		
        return std::move(level_buffers);
	}
	std::vector<std::vector<uint32_t>> get_level_buffer_dims(){
		return level_buffer_dims;
	}

private:
    bool use_sz = true;
	std::vector<std::vector<uint32_t>> level_dims;
	std::vector<std::vector<uint32_t>> level_buffer_dims;
	std::vector<std::vector<T>> level_buffers;
	std::vector<uint32_t> level_sizes;
	T * data_begin;
	long double sum=0;

	void init(const vector<size_t>& dims, size_t direction, size_t target_level){
		std::vector<uint32_t> dims_uint32;
		for(int i=0; i<dims.size(); i++){
			if(i != direction) dims_uint32.push_back(dims[i]);
		}
		level_dims = compute_level_dims_new(dims_uint32, target_level);
		level_sizes = compute_level_buffers_size_2D_coeff(level_dims, target_level, level_buffer_dims);
		for(int i=0; i<level_sizes.size(); i++){
			level_sizes[i] *= dims[direction];
			level_buffers.push_back(std::vector<T>(level_sizes[i]));
			// std::cout << "level_sizes[" << i << "] = " << level_sizes[i] << std::endl;
		}
	}

	size_t interleave_1D_level_0(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t count = 0;
		for(size_t i=begin; i<=end; i+=stride){
			// std::cout << (&(data[i]) - data_begin) << std::endl;
			buffer[count++] = data[i];
		}
		return count;
	}

	size_t compute_interpolant_difference_1D(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 1;
		size_t count = 0;
		for(i=1; i+1<n; i+=2){
			size_t c = begin + i * stride;
			sum += abs(data[c] - (data[c - stride] + data[c + stride]) / 2);
			// buffer[c] -= (buffer[c - stride] + buffer[c + stride]) / 2;
			// std::cout << (&data[c] - data_begin) << " -= (" << (&data[c - stride] - data_begin) << " + " << (&data[c + stride] - data_begin) << ")/2" << std::endl;
			buffer[count++] = data[c] - (data[c - stride] + data[c + stride]) / 2;
			data[c] = (data[c - stride] + data[c + stride]) / 2;
		}
		if(n % 2 == 0){
			size_t c = begin + (n - 1) * stride;
			sum += abs(data[c] - data[c - stride]);
			// std::cout << (&data[c] - data_begin) << " -= " << (&data[c - stride] - data_begin) << std::endl;
			// buffer[c] -= buffer[c - stride];
			buffer[count++] = data[c] - data[c - stride];
			data[c] = data[c - stride];
		}
		return count;
	}

	size_t compute_interpolant_difference_1D_diff_direct(const size_t begin, const size_t end, const size_t stride, const size_t interp_stride, bool last_even, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 0;
		size_t count = 0;
		if(!last_even){
			for(i=0; i<n; i+=2){
				size_t c = begin + i * stride;
				sum += abs(data[c] - (data[c - interp_stride] + data[c + interp_stride]) / 2);
				// std::cout << (&data[c] - data_begin) << " -= (" << (&data[c - interp_stride] - data_begin) << " + " << (&data[c + interp_stride] - data_begin) << ")/2" << std::endl;
				buffer[count++] = data[c] - (data[c - interp_stride] + data[c + interp_stride]) / 2;
				data[c] = (data[c - interp_stride] + data[c + interp_stride]) / 2;
			}
		} else {
			for(i=0; i<n; i+=2){
				size_t c = begin + i * stride;
				sum += abs(data[c] - data[c - interp_stride]);
				// std::cout << (&data[c] - data_begin) << " -= " << (&data[c - interp_stride] - data_begin) << std::endl;
				buffer[count++] = data[c] - data[c - interp_stride];
				data[c] = data[c - interp_stride];
			}
		}
		return count;
	}

	void decompose_interleave_level_1D_with_hierarchical_basis(T * data_pos, size_t n, T h, size_t current_level){
		size_t count_1D = 0;
		size_t buffer_index = current_level;
		count_1D = (current_level) ? compute_interpolant_difference_1D(0, n-1, h, data_pos, level_buffers[buffer_index].data()) : interleave_1D_level_0(0, n-1, h, data_pos, level_buffers[buffer_index].data());
		assert(count_1D == level_sizes[current_level]);
    }

	size_t interleave_2D_level_0(T * data_pos, size_t n1, size_t n2, size_t h){
		size_t count_2D = 0;
		size_t count;
		size_t stride_n2 = h;
		size_t stride_n1 = n2 * h;

		T * cur_data_pos = data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n2_begin = 0;
		size_t n2_end = n2 - 1;
		for(size_t i=0; i<n1; i+=h){
			count = interleave_1D_level_0(n2_begin, n2_end, stride_n2, cur_data_pos, cur_buffer_pos);
			cur_data_pos += stride_n1;
			cur_buffer_pos += count;
			count_2D += count;
		}
		return count_2D;
	}

	size_t compute_interpolant_difference_2D(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
		size_t count_2D = 0;
		size_t count;
		size_t stride_n1 = n2 * h;
		size_t stride_n2 = h;

		T * cur_data_pos;
		T * cur_buffer_pos;

		// old implementation: slow
		// cur_data_pos = data_pos;
		// cur_buffer_pos = level_buffers[((current_level - 1) * 2) + 1].data();
		// h = h << 1;
		// size_t n1_begin = 0;
		// size_t n1_end = (n1 - 1) * n2;
		// for(size_t i=0; i<n2; i+=h){
		// 	count = compute_interpolant_difference_1D(n1_begin, n1_end, stride_n1, cur_data_pos, cur_buffer_pos);
		// 	cur_data_pos += h;
		// 	cur_buffer_pos += count;
		// 	count_2D += count;
		// }

		// new implemenation: fast
		cur_data_pos = data_pos + stride_n1;
		cur_buffer_pos = level_buffers[((current_level - 1) * 2) + 1].data();
		// h = h << 1;
		size_t n1_begin = 0;
		size_t n1_end = (n1 - 1) * n2;
		size_t new_n1 = (n1_end - n1_begin) / stride_n1 + 1;
		size_t cur_begin = 0;
		size_t cur_end = n2 - 1;
		for(size_t i=1; i+1<new_n1; i+=2){
			count = compute_interpolant_difference_1D_diff_direct(cur_begin, cur_end, h, stride_n1, false, cur_data_pos, cur_buffer_pos);
			cur_data_pos += 2 * stride_n1;
			cur_buffer_pos += count;
			count_2D += count;
		}
		if(new_n1 % 2 == 0){
			count = compute_interpolant_difference_1D_diff_direct(cur_begin, cur_end, h, stride_n1, true, cur_data_pos, cur_buffer_pos);
			count_2D += count;
		}
		
		// h = h >> 1;
		cur_data_pos = data_pos;
		cur_buffer_pos = level_buffers[((current_level - 1) * 2) + 2].data();
		size_t n2_begin = 0;
		size_t n2_end = n2 - 1;
		for(size_t i=0; i<n1; i+=h){
			count = compute_interpolant_difference_1D(n2_begin, n2_end, stride_n2, cur_data_pos, cur_buffer_pos);
			cur_data_pos += n2 * h;
			cur_buffer_pos += count;
			count_2D += count;
		}
		
		return count_2D;
	}

	void decompose_interleave_2D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
        size_t count_2D = 0;
		count_2D = (current_level) ? compute_interpolant_difference_2D(data_pos, n1, n2, h, current_level) : interleave_2D_level_0(data_pos, n1, n2, h);
		assert(count_2D == (current_level) ? level_sizes[((current_level - 1) * 2) + 1] + level_sizes[((current_level - 1) * 2) + 2] : level_sizes[0]);
    }

	size_t interleave_3D_level_0(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n3 = h;
		size_t stride_n2 = n3*h;
		size_t stride_n1 = n2 * n3 * h;

		T * cur_data_pos = data_pos;
		T * temp_data_pos = cur_data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = interleave_1D_level_0(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3*h;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += n2*n3*h;
		}

		return count_3D;
	}

	size_t compute_interpolant_difference_3D(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n1 = n2 * n3 * h;
		size_t stride_n2 = n3 * h;
		size_t stride_n3 = h;

		T * cur_data_pos;
		T * temp_data_pos;
		T * cur_buffer_pos;
		size_t h2x = h << 1;

		size_t n1_begin = 0;
		size_t n1_end = (n1-1) * n2 * n3;
		size_t new_n1 = (n1_end - n1_begin) / stride_n1 + 1;
		size_t n2_begin = 0;
		size_t n2_end = (n2 - 1) * n3;
		size_t new_n2 = (n2_end - n2_begin) / stride_n2 + 1;
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		size_t new_n3 = (n3_end - n3_begin) / stride_n3 + 1;
		

		// old implementation: slow
		// cur_data_pos = data_pos;
		// cur_buffer_pos = level_buffers[((current_level - 1) * 3) + 1].data();
		// size_t n1_begin = 0;
		// size_t n1_end = (n1-1) * n2 * n3;
		// for(size_t i=0; i<n2; i+=h2x){
		// 	temp_data_pos = cur_data_pos;
		// 	for(size_t j=0; j<n3; j+=h2x){
		// 		count = compute_interpolant_difference_1D(n1_begin, n1_end, stride_n1, temp_data_pos, cur_buffer_pos);
		// 		temp_data_pos += h2x;
		// 		cur_buffer_pos += count;
		// 		count_3D += count;
		// 	}
		// 	cur_data_pos += n3 * h2x;
		// }

		// cur_data_pos = data_pos;
		// cur_buffer_pos = level_buffers[((current_level - 1) * 3) + 2].data();
		// size_t n2_begin = 0;
		// size_t n2_end = (n2-1) * n3;
		// for(size_t i=0; i<n1; i+=h){
		// 	temp_data_pos = cur_data_pos;
		// 	for(size_t j=0; j<n3; j+=h2x){
		// 		count = compute_interpolant_difference_1D(n2_begin, n2_end, stride_n2, temp_data_pos, cur_buffer_pos);
		// 		temp_data_pos += h2x;
		// 		cur_buffer_pos += count;
		// 		count_3D += count;
		// 	}
		// 	cur_data_pos += n2*n3*h;
		// }

		// new implementation: fast
		cur_data_pos = data_pos + stride_n1;
		cur_buffer_pos = level_buffers[((current_level - 1) * 3) + 1].data();
		for(size_t i=1; i+1<new_n1; i+=2){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h2x){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n1, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += 2 * stride_n1;
		}
		if(new_n1 % 2 == 0){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h2x){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n1, true, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
		}

		cur_data_pos = data_pos + stride_n2;
		cur_buffer_pos = level_buffers[((current_level - 1) * 3) + 2].data();
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=1; j+1<new_n2; j+=2){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n2, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			if(new_n2 % 2 == 0){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n2, true, temp_data_pos, cur_buffer_pos);
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += stride_n1;
		}
		

		cur_data_pos = data_pos;
		cur_buffer_pos = level_buffers[((current_level - 1) * 3) + 3].data();
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = compute_interpolant_difference_1D(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3*h;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += n2*n3*h;
		}
		return count_3D;
	}

	void decompose_interleave_3D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
        count_3D = (current_level != 0) ? compute_interpolant_difference_3D(data_pos, n1, n2, n3, h, current_level) : interleave_3D_level_0(data_pos, n1, n2, n3, h);
		assert(count_3D == (current_level) ? level_sizes[((current_level - 1) * 3) + 1] + level_sizes[((current_level - 1) * 3) + 2] + level_sizes[((current_level - 1) * 3) + 3] : level_sizes[0]);
    }

	size_t compute_interpolant_difference_1D_diff_direct_along_direction_2(const size_t begin, const size_t end, const size_t interp_stride, bool last_even, T * data, T * buffer){
		size_t n = (end - begin) + 1;
		size_t i = 0;
		size_t count = 0;
		if(!last_even){
			for(i=0; i<n; i++){
				size_t c = begin + i;
				sum += abs(data[c] - (data[c - interp_stride] + data[c + interp_stride]) / 2);
				// std::cout << (&data[c] - data_begin) << " -= (" << (&data[c - interp_stride] - data_begin) << " + " << (&data[c + interp_stride] - data_begin) << ")/2" << std::endl;
				// std::cout << data[c] - (data[c - interp_stride] + data[c + interp_stride]) / 2 << " = " << data[c] << " - (" << data[c - interp_stride] << " + " << data[c + interp_stride] << ")/2" << std::endl;
				buffer[count++] = data[c] - (data[c - interp_stride] + data[c + interp_stride]) / 2;
				data[c] = (data[c - interp_stride] + data[c + interp_stride]) / 2;
			}
		} else {
			for(i=0; i<n; i++){
				size_t c = begin + i;
				sum += abs(data[c] - data[c - interp_stride]);
				// std::cout << (&data[c] - data_begin) << " -= " << (&data[c - interp_stride] - data_begin) << std::endl;
				// std::cout << data[c] - data[c - interp_stride] << " = " << data[c] << " - " << data[c - interp_stride] << std::endl;
				buffer[count++] = data[c] - data[c - interp_stride];
				data[c] = data[c - interp_stride];
			}
		}
		return count;
	}

	size_t interleave_2D_level_0_along_direction_0(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n3 = h;
		size_t stride_n2 = n3*h;
		size_t stride_n1 = n2 * n3 * h;

		T * cur_data_pos = data_pos;
		T * temp_data_pos = cur_data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		for(size_t i=0; i<n1; i++){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = interleave_1D_level_0(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += n2*n3;
		}

		return count_3D;
	}

	size_t interleave_2D_level_0_along_direction_1(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n3 = h;
		size_t stride_n2 = n3*h;
		size_t stride_n1 = n2 * n3 * h;

		T * cur_data_pos = data_pos;
		T * temp_data_pos = cur_data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j++){
				count = interleave_1D_level_0(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += stride_n1;
		}

		return count_3D;
	}

	size_t interleave_2D_level_0_along_direction_2(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n3 = h;
		size_t stride_n2 = n3*h;
		size_t stride_n1 = n2 * n3 * h;

		T * cur_data_pos = data_pos;
		T * temp_data_pos = cur_data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = interleave_1D_level_0(n3_begin, n3_end, 1, temp_data_pos, cur_buffer_pos);
				temp_data_pos += stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += stride_n1;
		}

		return count_3D;
	}

	size_t compute_interpolant_difference_2D_along_direction_0(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n1 = n2 * n3 * h;
		size_t stride_n2 = n3 * h;
		size_t stride_n3 = h;

		T * cur_data_pos;
		T * temp_data_pos;
		T * cur_buffer_pos;
		size_t h2x = h << 1;

		size_t n1_begin = 0;
		size_t n1_end = (n1-1) * n2 * n3;
		size_t new_n1 = (n1_end - n1_begin) / stride_n1 + 1;
		size_t n2_begin = 0;
		size_t n2_end = (n2 - 1) * n3;
		size_t new_n2 = (n2_end - n2_begin) / stride_n2 + 1;
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		size_t new_n3 = (n3_end - n3_begin) / stride_n3 + 1;

		cur_data_pos = data_pos + stride_n2;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 1].data();
		for(size_t i=0; i<n1; i++){
			temp_data_pos = cur_data_pos;
			for(size_t j=1; j+1<new_n2; j+=2){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n2, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			if(new_n2 % 2 == 0){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n2, true, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += n2 * n3;
		}
		// std::cout << level_buffers[(current_level - 1) * 2 + 1].size() << " " << count_3D << std::endl;

		cur_data_pos = data_pos;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 2].data();
		for(size_t i=0; i<n1; i++){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = compute_interpolant_difference_1D(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += n2 * n3;
		}

		// std::cout << level_buffers[(current_level - 1) * 2 + 2].size() << " " << count_3D << std::endl;

		return count_3D;
	}

	size_t compute_interpolant_difference_2D_along_direction_1(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n1 = n2 * n3 * h;
		size_t stride_n2 = n3 * h;
		size_t stride_n3 = h;

		T * cur_data_pos;
		T * temp_data_pos;
		T * cur_buffer_pos;
		size_t h2x = h << 1;

		size_t n1_begin = 0;
		size_t n1_end = (n1-1) * n2 * n3;
		size_t new_n1 = (n1_end - n1_begin) / stride_n1 + 1;
		size_t n2_begin = 0;
		size_t n2_end = (n2 - 1) * n3;
		size_t new_n2 = (n2_end - n2_begin) / stride_n2 + 1;
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		size_t new_n3 = (n3_end - n3_begin) / stride_n3 + 1;

		cur_data_pos = data_pos + stride_n1;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 1].data();
		for(size_t i=1; i+1<new_n1; i+=2){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j++){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n1, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += 2 * stride_n1;
		}
		if(new_n1 % 2 == 0){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j++){
				count = compute_interpolant_difference_1D_diff_direct(n3_begin, n3_end, h, stride_n1, true, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3;
				cur_buffer_pos += count;
				count_3D += count;
			}
		}

		// std::cout << level_buffers[(current_level - 1) * 2 + 1].size() << " " << count_3D << std::endl;

		cur_data_pos = data_pos;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 2].data();
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j++){
				count = compute_interpolant_difference_1D(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
				temp_data_pos += n3;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += stride_n1;
		}

		// std::cout << level_buffers[(current_level - 1) * 2 + 2].size() << " " << count_3D << std::endl;

		return count_3D;
	}

	size_t compute_interpolant_difference_2D_along_direction_2(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
		size_t count;
		size_t stride_n1 = n2 * n3 * h;
		size_t stride_n2 = n3 * h;
		size_t stride_n3 = h;

		T * cur_data_pos;
		T * temp_data_pos;
		T * cur_buffer_pos;
		size_t h2x = h << 1;

		size_t n1_begin = 0;
		size_t n1_end = (n1-1) * n2 * n3;
		size_t new_n1 = (n1_end - n1_begin) / stride_n1 + 1;
		size_t n2_begin = 0;
		size_t n2_end = (n2 - 1) * n3;
		size_t new_n2 = (n2_end - n2_begin) / stride_n2 + 1;
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		size_t new_n3 = (n3_end - n3_begin) / stride_n3 + 1;

		cur_data_pos = data_pos + stride_n1;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 1].data();
		for(size_t i=1; i+1<new_n1; i+=2){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h2x){
				count = compute_interpolant_difference_1D_diff_direct_along_direction_2(n3_begin, n3_end, stride_n1, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += 2 * stride_n1;
		}
		if(new_n1 % 2 == 0){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h2x){
				count = compute_interpolant_difference_1D_diff_direct_along_direction_2(n3_begin, n3_end, stride_n1, true, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
		}

		// std::cout << level_buffers[(current_level - 1) * 2 + 1].size() << " " << count_3D << std::endl;

		cur_data_pos = data_pos + stride_n2;
		cur_buffer_pos = level_buffers[(current_level - 1) * 2 + 2].data();
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=1; j+1<new_n2; j+=2){
				count = compute_interpolant_difference_1D_diff_direct_along_direction_2(n3_begin, n3_end, stride_n2, false, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			if(new_n2 % 2 == 0){
				count = compute_interpolant_difference_1D_diff_direct_along_direction_2(n3_begin, n3_end, stride_n2, true, temp_data_pos, cur_buffer_pos);
				temp_data_pos += 2 * stride_n2;
				cur_buffer_pos += count;
				count_3D += count;
			}
			cur_data_pos += stride_n1;
		}

		// std::cout << level_buffers[(current_level - 1) * 2 + 2].size() << " " << count_3D << std::endl;

		return count_3D;
	}

	void decompose_interleave_3D_HB_with_direction(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t direction, size_t current_level){
		size_t count_3D = 0;
		switch (direction){
			case 0:
			{
				if(current_level == 0) count_3D = interleave_2D_level_0_along_direction_0(data_pos, n1, n2, n3, h);
				else count_3D = compute_interpolant_difference_2D_along_direction_0(data_pos, n1, n2, n3, h, current_level);
				break;
			}
			case 1:
			{
				if(current_level == 0) count_3D = interleave_2D_level_0_along_direction_1(data_pos, n1, n2, n3, h);
				else count_3D = compute_interpolant_difference_2D_along_direction_1(data_pos, n1, n2, n3, h, current_level);
				break;
			}
			case 2:
			{
				if(current_level == 0) count_3D = interleave_2D_level_0_along_direction_2(data_pos, n1, n2, n3, h);
				else count_3D = compute_interpolant_difference_2D_along_direction_2(data_pos, n1, n2, n3, h, current_level);
				break;
			}
			default:
				std::cerr << "Unsupported direction" << std::endl;
				exit(-1);
		}
		// std::cout << "current_level = " << current_level << std::endl;
		// std::cout << "count_3D = " << count_3D << ", level_size = " << ((current_level) ? level_sizes[((current_level - 1) * 2) + 1] + level_sizes[((current_level - 1) * 2) + 2] : level_sizes[0]) << std::endl;
		assert(count_3D == (current_level) ? level_sizes[((current_level - 1) * 2) + 1] + level_sizes[((current_level - 1) * 2) + 2] : level_sizes[0]);
	}
};
}

#endif