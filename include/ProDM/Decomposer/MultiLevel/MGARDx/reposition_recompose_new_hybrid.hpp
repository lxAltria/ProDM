#ifndef _MGARD_REPOSITION_RECOMPOSE_NEW_HYBRID_HPP
#define _MGARD_REPOSITION_RECOMPOSE_NEW_HYBRID_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "reorder.hpp"
#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "correction.hpp"

namespace MGARD{

using namespace std;

// Hybrid interpolation: current_level == 1 uses linear, current_level >= 2 uses cubic.
// current_level == 0 always repositions coarse grid data (unchanged).
template <class T>
class Repositioner_Recomposer_new_hybrid{
public:
	Repositioner_Recomposer_new_hybrid(bool use_sz_=true){
            use_sz = use_sz_;
        };
	~Repositioner_Recomposer_new_hybrid(){
	};
	std::vector<T> recompose(std::vector<std::vector<T>>& level_buffers_, const vector<size_t>& dims, size_t target_level, vector<size_t> strides=vector<size_t>()){
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}

		init(dims, target_level);
		data_buffer.resize(num_elements);
		std::fill(data_buffer.begin(), data_buffer.end(), 0);
		data_begin = data_buffer.data();

		level_buffers = level_buffers_;
		if(dims.size() == 1){
			size_t h = 1 << target_level;
			size_t n = dims[0];
			for(int current_level=0; current_level <= target_level; current_level++){
				bool use_linear = (current_level == 1);
				if(use_linear) recompose_level_1D_with_hierarchical_basis(data_buffer.data(), n, h, current_level);
				else           recompose_level_1D_cubic_with_hierarchical_basis(data_buffer.data(), n, h, current_level);
				h >>= 1;
			}
		}
		else if(dims.size() == 2){
			size_t h = 1 << target_level;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			for(int current_level=0; current_level <= target_level; current_level++){
				bool use_linear = (current_level == 1);
				if(use_linear) recompose_level_2D_with_hierarchical_basis(data_buffer.data(), n1, n2, h, current_level);
				else           recompose_level_2D_cubic_with_hierarchical_basis(data_buffer.data(), n1, n2, h, current_level);
				h >>= 1;
			}
		}
		else if(dims.size() == 3){
			size_t h = 1 << target_level;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			size_t n3 = dims[2];
			for(int current_level=0; current_level <= target_level; current_level++){
				bool use_linear = (current_level == 1);
				if(use_linear) recompose_level_3D_with_hierarchical_basis(data_buffer.data(), n1, n2, n3, h, current_level);
				else           recompose_level_3D_cubic_with_hierarchical_basis(data_buffer.data(), n1, n2, n3, h, current_level);
				h >>= 1;
			}
		}
        return std::move(data_buffer);
	}
	std::vector<std::vector<uint32_t>> get_level_buffer_dims(){
		return level_buffer_dims;
	}
	std::vector<uint32_t> interp_order = {0, 1, 2};

private:
    bool use_sz = true;
	std::vector<std::vector<T>> level_buffers;
	std::vector<std::vector<uint32_t>> level_buffer_dims;
	std::vector<uint32_t> level_sizes;
	std::vector<T> data_buffer;
	T * data_begin;

	void init(const vector<size_t>& dims, size_t target_level){
		std::vector<uint32_t> dims_uint32(dims.begin(), dims.end());
		auto level_dims = compute_level_dims_new(dims_uint32, target_level);
		level_sizes = compute_level_buffers_size_generic(level_dims, target_level, interp_order, level_buffer_dims);
	}

	size_t reposition_1D_level_0(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t count = 0;
		for(size_t i=begin; i<=end; i+=stride){
			data[i] = buffer[count++];
		}
		return count;
	}

	size_t compute_interpolant_difference_1D(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 1;
		size_t count = 0;
		for(i=1; i+1<n; i+=2){
			size_t c = begin + i * stride;
			data[c] = (data[c - stride] + data[c + stride]) / 2 + buffer[count];
			count++;
		}
		if(n % 2 == 0){
			size_t c = begin + (n - 1) * stride;
			data[c] = data[c - stride] + buffer[count];
			count++;
		}
		return count;
	}

	size_t compute_interpolant_difference_1D_cubic(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 1;
		size_t count = 0;
		size_t c;
		size_t stride3x = 3 * stride;
		size_t stride5x = 5 * stride;
		c = begin + i * stride;
		data[c] = interp_quad_1(data[c - stride], data[c + stride], data[c + stride3x]) + buffer[count];
		count++;
		for(i=3; i+3<n; i+=2){
			c = begin + i * stride;
			data[c] = interp_cubic(data[c - stride3x], data[c - stride], data[c + stride], data[c + stride3x]) + buffer[count];
			count++;
		}
		c = begin + i * stride;
		data[c] = interp_quad_2(data[c - stride3x], data[c - stride], data[c + stride]) + buffer[count];
		count++;
		if(!(n&1)){
			c = begin + (n - 1) * stride;
			data[c] = interp_quad_3(data[c - stride5x], data[c - stride3x], data[c - stride]) + buffer[count];
			count++;
		}
		return count;
	}

	size_t compute_interpolant_difference_1D_diff_direct(const size_t begin, const size_t end, const size_t stride, const size_t interp_stride, bool last_even, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 0;
		size_t count = 0;
		if(!last_even){
			for(i=0; i<n; i++){
				size_t c = begin + i * stride;
				data[c] = (data[c - interp_stride] + data[c + interp_stride]) / 2 + buffer[count];
				count++;
			}
		} else {
			for(i=0; i<n; i++){
				size_t c = begin + i * stride;
				data[c] = data[c - interp_stride] + buffer[count];
				count++;
			}
		}
		return count;
	}

	size_t compute_interpolant_difference_1D_cubic_diff_direct(const size_t begin, const size_t end, const size_t stride, const size_t interp_stride, int interp_type, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 0;
		size_t count = 0;
		size_t c;
		size_t interp_stride3x = 3 * interp_stride;
		size_t interp_stride5x = 5 * interp_stride;
		switch (interp_type)
		{
			case FIRST_LINE:
			{
				for(i=0; i<n; i++){
					c = begin + i * stride;
					data[c] = interp_quad_1(data[c - interp_stride], data[c + interp_stride], data[c + interp_stride3x]) + buffer[count];
					count++;
				}
				break;
			}
			case NORMAL_LINE:
			{
				for(i=0; i<n; i++){
					c = begin + i * stride;
					data[c] = interp_cubic(data[c - interp_stride3x], data[c - interp_stride], data[c + interp_stride], data[c + interp_stride3x]) + buffer[count];
					count++;
				}
				break;
			}
			case LAST_LINE:
			{
				for(i=0; i<n; i++){
					c = begin + i * stride;
					data[c] = interp_quad_2(data[c - interp_stride3x], data[c - interp_stride], data[c + interp_stride]) + buffer[count];
					count++;
				}
				break;
			}
			case LAST_EVEN_LINE:
			{
				for(i=0; i<n; i++){
					c = begin + i * stride;
					data[c] = interp_quad_3(data[c - interp_stride5x], data[c - interp_stride3x], data[c - interp_stride]) + buffer[count];
					count++;
				}
				break;
			}
			default:
				std::cout << "No cubic interpolation is suitable" << std::endl;
				exit(-1);
				break;
		}
		return count;
	}

	size_t recover_from_interpolant_difference_1D(const size_t begin, const size_t end, const size_t stride, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 1;
		size_t count = 0;
		for(i=1; i+1<n; i+=2){
			size_t c = begin + i * stride;
			data[c] += buffer[count++];
		}
		if(n % 2 == 0){
			size_t c = begin + (n - 1) * stride;
			data[c] += buffer[count++];
		}
		return count;
	}

	size_t recover_from_interpolant_difference_1D_diff_direct(const size_t begin, const size_t end, const size_t stride, const size_t interp_stride, bool last_even, T * data, T * buffer){
		size_t n = (end - begin) / stride + 1;
		size_t i = 0;
		size_t count = 0;
		if(!last_even){
			for(i=0; i<n; i++){
				size_t c = begin + i * stride;
				data[c] += buffer[count++];
			}
		} else {
			for(i=0; i<n; i++){
				size_t c = begin + i * stride;
				data[c] += buffer[count++];
			}
		}
		return count;
	}

	void recompose_level_1D_with_hierarchical_basis(T * data_pos, size_t n, T h, size_t current_level){
		size_t count_1D = 0;
		if(current_level){
			count_1D = compute_interpolant_difference_1D(0, n-1, h, data_pos, level_buffers[current_level].data());
		}
		else{
			count_1D = reposition_1D_level_0(0, n-1, h, data_pos, level_buffers[current_level].data());
		}
		assert(count_1D == level_sizes[current_level]);
    }

	void recompose_level_1D_cubic_with_hierarchical_basis(T * data_pos, size_t n, T h, size_t current_level){
		size_t count_1D = 0;
		if(current_level){
			count_1D = compute_interpolant_difference_1D_cubic(0, n-1, h, data_pos, level_buffers[current_level].data());
		}
		else{
			count_1D = reposition_1D_level_0(0, n-1, h, data_pos, level_buffers[current_level].data());
		}
		assert(count_1D == level_sizes[current_level]);
    }

	size_t reposition_2D_level_0(T * data_pos, size_t n1, size_t n2, size_t h){
		size_t count_2D = 0;
		size_t count;
		size_t stride_n2 = h;
		size_t stride_n1 = n2 * h;

		T * cur_data_pos = data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n2_begin = 0;
		size_t n2_end = n2 - 1;
		for(size_t i=0; i<n1; i+=h){
			count = reposition_1D_level_0(n2_begin, n2_end, stride_n2, cur_data_pos, cur_buffer_pos);
			cur_data_pos += stride_n1;
			cur_buffer_pos += count;
			count_2D += count;
		}
		return count_2D;
	}

	size_t compute_interpolant_difference_2D(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
		size_t count_2D = 0;
		size_t count;
		size_t h2x = h << 1;
		size_t stride[2] = {n2*h, h};
		size_t new_n[2] = {(n1-1)/h + 1, (n2-1)/h + 1};
		bool interpolated[2] = {false, false};
		T * d0_pos;
		T * cur_buffer_pos;

		for(int s = 0; s < 2; s++){
			uint32_t id = interp_order[s];
			cur_buffer_pos = level_buffers[((current_level - 1) * 2) + s + 1].data();
			size_t step0 = interpolated[0] ? h : h2x;
			size_t dstep0 = n2 * step0;
			size_t y_stride = interpolated[1] ? h : h2x;

			if(id == 1){
				d0_pos = data_pos;
				for(size_t i=0; i<n1; i+=step0){
					count = compute_interpolant_difference_1D(0, n2-1, h, d0_pos, cur_buffer_pos);
					d0_pos += dstep0;
					cur_buffer_pos += count;
					count_2D += count;
				}
			} else {
				d0_pos = data_pos + stride[0];
				for(size_t i=1; i+1<new_n[0]; i+=2){
					count = compute_interpolant_difference_1D_diff_direct(0, n2-1, y_stride, stride[0], false, d0_pos, cur_buffer_pos);
					d0_pos += 2 * stride[0];
					cur_buffer_pos += count;
					count_2D += count;
				}
				if(!(new_n[0] & 1)){
					count = compute_interpolant_difference_1D_diff_direct(0, n2-1, y_stride, stride[0], true, d0_pos, cur_buffer_pos);
					count_2D += count;
				}
			}
			interpolated[id] = true;
		}
		return count_2D;
	}

	size_t compute_interpolant_difference_2D_cubic(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
		size_t count_2D = 0;
		size_t count;
		size_t h2x = h << 1;
		size_t stride[2] = {n2*h, h};
		size_t new_n[2] = {(n1-1)/h + 1, (n2-1)/h + 1};
		bool interpolated[2] = {false, false};
		T * d0_pos;
		T * cur_buffer_pos;

		for(int s = 0; s < 2; s++){
			uint32_t id = interp_order[s];
			cur_buffer_pos = level_buffers[((current_level - 1) * 2) + s + 1].data();
			size_t step0 = interpolated[0] ? h : h2x;
			size_t dstep0 = n2 * step0;
			size_t y_stride = interpolated[1] ? h : h2x;

			if(id == 1){
				d0_pos = data_pos;
				for(size_t i=0; i<n1; i+=step0){
					count = compute_interpolant_difference_1D_cubic(0, n2-1, h, d0_pos, cur_buffer_pos);
					d0_pos += dstep0;
					cur_buffer_pos += count;
					count_2D += count;
				}
			} else {
				d0_pos = data_pos + stride[0];
				// First
				count = compute_interpolant_difference_1D_cubic_diff_direct(0, n2-1, y_stride, stride[0], FIRST_LINE, d0_pos, cur_buffer_pos);
				d0_pos += 2 * stride[0];
				cur_buffer_pos += count;
				count_2D += count;
				// Normal
				for(size_t i=3; i+3<new_n[0]; i+=2){
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n2-1, y_stride, stride[0], NORMAL_LINE, d0_pos, cur_buffer_pos);
					d0_pos += 2 * stride[0];
					cur_buffer_pos += count;
					count_2D += count;
				}
				// Last
				count = compute_interpolant_difference_1D_cubic_diff_direct(0, n2-1, y_stride, stride[0], LAST_LINE, d0_pos, cur_buffer_pos);
				d0_pos += 2 * stride[0];
				cur_buffer_pos += count;
				count_2D += count;
				// Last even
				if(!(new_n[0] & 1)){
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n2-1, y_stride, stride[0], LAST_EVEN_LINE, d0_pos, cur_buffer_pos);
					count_2D += count;
				}
			}
			interpolated[id] = true;
		}
		return count_2D;
	}

	void recompose_level_2D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
        size_t count_2D = 0;
		if(current_level){
			count_2D = compute_interpolant_difference_2D(data_pos, n1, n2, h, current_level);
		}
		else{
			count_2D = reposition_2D_level_0(data_pos, n1, n2, h);
		}
		assert(count_2D == (current_level) ? level_sizes[((current_level - 1) * 2) + 1] + level_sizes[((current_level - 1) * 2) + 2] : level_sizes[0]);
    }

	void recompose_level_2D_cubic_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t h, size_t current_level){
        size_t count_2D = 0;
		if(current_level){
			count_2D = compute_interpolant_difference_2D_cubic(data_pos, n1, n2, h, current_level);
		}
		else{
			count_2D = reposition_2D_level_0(data_pos, n1, n2, h);
		}
		assert(count_2D == (current_level) ? level_sizes[((current_level - 1) * 2) + 1] + level_sizes[((current_level - 1) * 2) + 2] : level_sizes[0]);
    }

	size_t reposition_3D_level_0(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h){
		size_t count_3D = 0;
		size_t count;
		T * cur_data_pos = data_pos;
		T * temp_data_pos = cur_data_pos;
		T * cur_buffer_pos = level_buffers[0].data();
		size_t n3_begin = 0;
		size_t n3_end = n3 - 1;
		for(size_t i=0; i<n1; i+=h){
			temp_data_pos = cur_data_pos;
			for(size_t j=0; j<n2; j+=h){
				count = reposition_1D_level_0(n3_begin, n3_end, h, temp_data_pos, cur_buffer_pos);
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
		size_t h2x = h << 1;
		size_t stride[3] = {n2*n3*h, n3*h, h};
		size_t new_n[3] = {(n1-1)/h + 1, (n2-1)/h + 1, (n3-1)/h + 1};
		bool interpolated[3] = {false, false, false};
		T * d0_pos;
		T * d1_pos;
		T * cur_buffer_pos;

		for(int s = 0; s < 3; s++){
			uint32_t id = interp_order[s];
			cur_buffer_pos = level_buffers[((current_level - 1) * 3) + s + 1].data();
			size_t step0 = interpolated[0] ? h : h2x;
			size_t step1 = interpolated[1] ? h : h2x;
			size_t dstep0 = n2 * n3 * step0;
			size_t dstep1 = n3 * step1;
			size_t z_stride = interpolated[2] ? h : h2x;

			if(id == 2){
				d0_pos = data_pos;
				for(size_t i=0; i<n1; i+=step0){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D(0, n3-1, h, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += dstep0;
				}
			} else if(id == 0){
				d0_pos = data_pos + stride[0];
				for(size_t i=1; i+1<new_n[0]; i+=2){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D_diff_direct(0, n3-1, z_stride, stride[0], false, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += 2 * stride[0];
				}
				if(!(new_n[0] & 1)){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D_diff_direct(0, n3-1, z_stride, stride[0], true, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
				}
			} else {
				d0_pos = data_pos + stride[1];
				for(size_t i=0; i<n1; i+=step0){
					d1_pos = d0_pos;
					for(size_t j=1; j+1<new_n[1]; j+=2){
						count = compute_interpolant_difference_1D_diff_direct(0, n3-1, z_stride, stride[1], false, d1_pos, cur_buffer_pos);
						d1_pos += 2 * stride[1];
						cur_buffer_pos += count;
						count_3D += count;
					}
					if(!(new_n[1] & 1)){
						count = compute_interpolant_difference_1D_diff_direct(0, n3-1, z_stride, stride[1], true, d1_pos, cur_buffer_pos);
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += dstep0;
				}
			}
			interpolated[id] = true;
		}
		return count_3D;
	}

	size_t compute_interpolant_difference_3D_cubic(T * data_pos, size_t n1, size_t n2, size_t n3, size_t h, size_t current_level){
		size_t count_3D = 0;
		size_t count;
		size_t h2x = h << 1;
		size_t stride[3] = {n2*n3*h, n3*h, h};
		size_t new_n[3] = {(n1-1)/h + 1, (n2-1)/h + 1, (n3-1)/h + 1};
		bool interpolated[3] = {false, false, false};
		T * d0_pos;
		T * d1_pos;
		T * cur_buffer_pos;

		for(int s = 0; s < 3; s++){
			uint32_t id = interp_order[s];
			cur_buffer_pos = level_buffers[((current_level - 1) * 3) + s + 1].data();
			size_t step0 = interpolated[0] ? h : h2x;
			size_t step1 = interpolated[1] ? h : h2x;
			size_t dstep0 = n2 * n3 * step0;
			size_t dstep1 = n3 * step1;
			size_t z_stride = interpolated[2] ? h : h2x;

			if(id == 2){
				d0_pos = data_pos;
				for(size_t i=0; i<n1; i+=step0){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D_cubic(0, n3-1, h, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += dstep0;
				}
			} else if(id == 0){
				d0_pos = data_pos + stride[0];
				// First plane
				d1_pos = d0_pos;
				for(size_t j=0; j<n2; j+=step1){
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[0], FIRST_LINE, d1_pos, cur_buffer_pos);
					d1_pos += dstep1;
					cur_buffer_pos += count;
					count_3D += count;
				}
				d0_pos += 2 * stride[0];
				// Normal planes
				for(size_t i=3; i+3<new_n[0]; i+=2){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[0], NORMAL_LINE, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += 2 * stride[0];
				}
				// Last plane
				d1_pos = d0_pos;
				for(size_t j=0; j<n2; j+=step1){
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[0], LAST_LINE, d1_pos, cur_buffer_pos);
					d1_pos += dstep1;
					cur_buffer_pos += count;
					count_3D += count;
				}
				d0_pos += 2 * stride[0];
				// Last even plane
				if(!(new_n[0] & 1)){
					d1_pos = d0_pos;
					for(size_t j=0; j<n2; j+=step1){
						count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[0], LAST_EVEN_LINE, d1_pos, cur_buffer_pos);
						d1_pos += dstep1;
						cur_buffer_pos += count;
						count_3D += count;
					}
				}
			} else {
				d0_pos = data_pos + stride[1];
				for(size_t i=0; i<n1; i+=step0){
					d1_pos = d0_pos;
					// First line
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[1], FIRST_LINE, d1_pos, cur_buffer_pos);
					d1_pos += 2 * stride[1];
					cur_buffer_pos += count;
					count_3D += count;
					// Normal lines
					for(size_t j=3; j+3<new_n[1]; j+=2){
						count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[1], NORMAL_LINE, d1_pos, cur_buffer_pos);
						d1_pos += 2 * stride[1];
						cur_buffer_pos += count;
						count_3D += count;
					}
					// Last line
					count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[1], LAST_LINE, d1_pos, cur_buffer_pos);
					d1_pos += 2 * stride[1];
					cur_buffer_pos += count;
					count_3D += count;
					// Last even line
					if(!(new_n[1] & 1)){
						count = compute_interpolant_difference_1D_cubic_diff_direct(0, n3-1, z_stride, stride[1], LAST_EVEN_LINE, d1_pos, cur_buffer_pos);
						cur_buffer_pos += count;
						count_3D += count;
					}
					d0_pos += dstep0;
				}
			}
			interpolated[id] = true;
		}
		return count_3D;
	}

	void recompose_level_3D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t current_level){
		size_t count_3D = 0;
		if(current_level){
			count_3D = compute_interpolant_difference_3D(data_pos, n1, n2, n3, h, current_level);
		}
		else{
			count_3D = reposition_3D_level_0(data_pos, n1, n2, n3, h);
		}
		assert(count_3D == (current_level) ? level_sizes[((current_level - 1) * 3) + 1] + level_sizes[((current_level - 1) * 3) + 2] + level_sizes[((current_level - 1) * 3) + 3] : level_sizes[0]);
    }

	void recompose_level_3D_cubic_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t current_level){
		size_t count_3D = 0;
		if(current_level){
			count_3D = compute_interpolant_difference_3D_cubic(data_pos, n1, n2, n3, h, current_level);
		}
		else{
			count_3D = reposition_3D_level_0(data_pos, n1, n2, n3, h);
		}
		assert(count_3D == (current_level) ? level_sizes[((current_level - 1) * 3) + 1] + level_sizes[((current_level - 1) * 3) + 2] + level_sizes[((current_level - 1) * 3) + 3] : level_sizes[0]);
    }
};
}

#endif
