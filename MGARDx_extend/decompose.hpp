#ifndef _MGARD_DECOMPOSE_HPP
#define _MGARD_DECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "reorder.hpp"
#include "utils.hpp"
#include "correction.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Decomposer{
public:
	Decomposer(bool use_sz_=true){
            use_sz = use_sz_;
        };
	~Decomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
    // return levels
	int decompose(T * data_, const vector<size_t>& dims, size_t target_level, bool hierarchical=false, bool cubic=false, vector<size_t> strides=vector<size_t>()){
		data = data_;
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
		if(strides.size() == 0){
			strides = vector<size_t>(dims.size());
			size_t stride = 1;
			for(int i=dims.size()-1; i>=0; i--){
				strides[i] = stride;
				stride *= dims[i];
			}
			for(int i=0; i<dims.size(); i++){
				cout << strides[i] << " ";
			}
			cout << endl;
		}
		data_buffer_size = num_elements * sizeof(T);
        int max_level = log2(*min_element(dims.begin(), dims.end()));
        if(target_level > max_level) target_level = max_level;
		init(dims);
        current_dims.resize(dims.size());
		if(dims.size() == 1){
			size_t h = 1;
			size_t n = dims[0];
			for(int i=0; i<target_level; i++){
				if (cubic)
					decompose_level_1D_cubic(data, n, h);
				else
					hierarchical ? decompose_level_1D_with_hierarchical_basis(data, n, h) : decompose_level_1D(data, n, h);
				n = (n >> 1) + 1;
				h <<= 1;
			}
		}
		else if(dims.size() == 2){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			for(int i=0; i<target_level; i++){
				if (cubic)
					decompose_level_2D_cubic(data, n1, n2, (T)h, strides[0]);
				else
					hierarchical ? decompose_level_2D_with_hierarchical_basis(data, n1, n2, (T)h, strides[0]) : decompose_level_2D(data, n1, n2, (T)h, strides[0]);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				h <<= 1;
			}
		}
		else if(dims.size() == 3){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			size_t n3 = dims[2];
			for(int i=0; i<target_level; i++){
                current_dims[0] = (n1 >> 1) + 1;
                current_dims[1] = (n2 >> 1) + 1;
                current_dims[2] = (n3 >> 1) + 1;
				if (cubic)
					decompose_level_3D_cubic(data, n1, n2, n3, (T)h, strides[0], strides[1]);
				else
					hierarchical ? decompose_level_3D_with_hierarchical_basis(data, n1, n2, n3, (T)h, strides[0], strides[1]) : decompose_level_3D(data, n1, n2, n3, (T)h, strides[0], strides[1]);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				n3 = (n3 >> 1) + 1;
				h <<= 1;
			}
		}
        return target_level;
	}

private:
	unsigned int default_batch_size = 32;
	size_t data_buffer_size = 0;
    bool use_sz = true;
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;
    vector<size_t> current_dims;

	void init(const vector<size_t>& dims){
		size_t buffer_size = default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		// cerr << "buffer_size = " << buffer_size << endl;
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(data_buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l
	// overwrite the data in N_l \ N_(l-1) in place
	void compute_interpolant_difference_1D(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] -= (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}
	void compute_interpolant_difference_1D_cubic(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
		coeff_buffer[0] -= interp_quad_1(nodal_buffer[0], nodal_buffer[1], nodal_buffer[2]);
		coeff_buffer[n_coeff-1] -= interp_quad_2(nodal_buffer[n_coeff-2], nodal_buffer[n_coeff-1], nodal_buffer[n_coeff]);
		for(int i=1; i+1<n_coeff; i++){
			coeff_buffer[i] -= interp_cubic(nodal_buffer[i-1], nodal_buffer[i], nodal_buffer[i+1], nodal_buffer[i+2]); 
		}
	}

	void add_correction(size_t n_nodal, T * nodal_buffer){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] += correction_buffer[i];
		}
	}
	// decompose a level with n element and the given stride
	// to a level with n/2 element
	void decompose_level_1D(T * data_pos, size_t n, T h, bool nodal_row=true){
		size_t n_nodal = (n >> 1) + 1;
		size_t n_coeff = n - n_nodal;
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
		compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
		if(nodal_row) compute_load_vector_nodal_row(load_v_buffer, n_nodal, n_coeff, h, coeff_buffer);
        else compute_load_vector_coeff_row(load_v_buffer, n_nodal, n_coeff, h, nodal_buffer, coeff_buffer);
		compute_correction(correction_buffer, n_nodal, h, load_v_buffer);
		add_correction(n_nodal, nodal_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}
    void decompose_level_1D_with_hierarchical_basis(T * data_pos, size_t n, T h, bool nodal_row=true){
        size_t n_nodal = (n >> 1) + 1;
        size_t n_coeff = n - n_nodal;
        T * nodal_buffer = data_buffer;
        T * coeff_buffer = data_buffer + n_nodal;
        data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
        compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
    }
	void decompose_level_1D_cubic(T * data_pos, size_t n, T h, bool nodal_row=true){
		size_t n_nodal = (n >> 1) + 1;
        size_t n_coeff = n - n_nodal;
        T * nodal_buffer = data_buffer;
        T * coeff_buffer = data_buffer + n_nodal;
        data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
        compute_interpolant_difference_1D_cubic(n_coeff, nodal_buffer, coeff_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l for the coefficient rows in 2D
	// overwrite the data in N_l \ N_(l-1) in place
	// Note: interpolant difference in the nodal rows have already been computed
	void compute_interpolant_difference_2D_vertical(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		bool even_n2 = !(n2 & 1);
		T * n1_nodal_data = data_pos;
		T * n1_coeff_data = data_pos + n1_nodal * stride;
		for(int i=0; i<n1_coeff; i++){
            const T * nodal_pos = n1_nodal_data + i * stride;
            T * coeff_pos = n1_coeff_data + i * stride;
            // TODO: optimize average computation
            T * nodal_coeff_pos = coeff_pos;	// coeffcients in nodal rows
            T * coeff_coeff_pos = coeff_pos + n2_nodal;	// coefficients in coeffcients rows
            for(int j=0; j<n2_coeff; j++){
                // coefficients in nodal columns
                *(nodal_coeff_pos++) -= (nodal_pos[j] + nodal_pos[stride + j]) / 2;
                // coefficients in centers
                *(coeff_coeff_pos++) -= (nodal_pos[j] + nodal_pos[j + 1] + nodal_pos[stride + j] + nodal_pos[stride + j + 1]) / 4;
            }
            // compute the last (or second last if n2 is even) nodal column
            *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff] + nodal_pos[stride + n2_coeff]) / 2;
            if(even_n2){
                // compute the last nodal column
                *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff + 1] + nodal_pos[stride + n2_coeff + 1]) / 2;
            }
		}
	}
	void compute_interpolant_difference_2D_vertical_cubic(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		bool even_n2 = !(n2 & 1);
		T * n1_nodal_data = data_pos;
		T * n1_coeff_data = data_pos + n1_nodal * stride;
		for(int i=0; i<n1_coeff; i++){
            const T * nodal_pos = n1_nodal_data + i * stride;
            T * coeff_pos = n1_coeff_data + i * stride;	
            // TODO: optimize average computation
            T * nodal_coeff_pos = coeff_pos;	// coeffcients in nodal rows
            T * coeff_coeff_pos = coeff_pos + n2_nodal;	// coefficients in coeffcients rows
			if (i == 0){	// when i = 0 which means the first coeff row
				// when j = 0 which means the first coeff element in the first coeff row
				*(nodal_coeff_pos++) -= interp_quad_1(nodal_pos[0], nodal_pos[stride], nodal_pos[2*stride]);
				*(coeff_coeff_pos++) -= interp_quad_1(interp_quad_1(nodal_pos[0], nodal_pos[stride], nodal_pos[2*stride]), 
														interp_quad_1(nodal_pos[1], nodal_pos[stride + 1], nodal_pos[2*stride + 1]),
															interp_quad_1(nodal_pos[2], nodal_pos[stride + 2], nodal_pos[2*stride + 2]));
				for(int j=1; j+1<n2_coeff; j++){	// regular case
					*(nodal_coeff_pos++) -= interp_quad_1(nodal_pos[j], nodal_pos[stride + j], nodal_pos[2 * stride + j]);
					*(coeff_coeff_pos++) -= interp_cubic(interp_quad_1(nodal_pos[j - 1], nodal_pos[stride + j - 1], nodal_pos[2*stride + j - 1]),
															interp_quad_1(nodal_pos[j], nodal_pos[stride + j], nodal_pos[2*stride + j]),
																interp_quad_1(nodal_pos[j + 1], nodal_pos[stride + j + 1], nodal_pos[2*stride + j + 1]),
																	interp_quad_1(nodal_pos[j + 2], nodal_pos[stride + j + 2], nodal_pos[2*stride + j + 2]));
				}
				// when j = n2_coeff - 1 which means the last coeff element in the first coeff row
				*(nodal_coeff_pos++) -= interp_quad_1(nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1], nodal_pos[2*stride + n2_coeff - 1]);
				*(coeff_coeff_pos++) -= interp_quad_2(interp_quad_1(nodal_pos[n2_coeff - 2], nodal_pos[stride + n2_coeff - 2], nodal_pos[2*stride + n2_coeff - 2]),
														interp_quad_1(nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1], nodal_pos[2*stride + n2_coeff - 1]),
															interp_quad_1(nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff], nodal_pos[2*stride + n2_coeff]));
				// compute the last (or second last if n2 is even) nodal column
				*(nodal_coeff_pos++) -= interp_quad_1(nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff], nodal_pos[2 * stride + n2_coeff]);
				if (even_n2){	
					// compute the last nodal column
					*(nodal_coeff_pos++) -= interp_quad_1(nodal_pos[n2_coeff + 1], nodal_pos[stride + n2_coeff + 1], nodal_pos[2*stride + n2_coeff + 1]);
				}
			}
			else if (i+1 == n1_coeff){	// when i+1 = n1_coeff which means the last coeff row
				// when j = 0 which means the first coeff element in the first coeff row
				*(nodal_coeff_pos++) -= interp_quad_2(nodal_pos[-stride], nodal_pos[0], nodal_pos[stride]);
				*(coeff_coeff_pos++) -= interp_quad_1(interp_quad_2(nodal_pos[-stride], nodal_pos[0], nodal_pos[stride]),
														interp_quad_2(nodal_pos[-stride + 1], nodal_pos[1], nodal_pos[stride + 1]),
															interp_quad_2(nodal_pos[-stride + 2], nodal_pos[2], nodal_pos[stride + 2]));
				for(int j=1; j+1<n2_coeff; j++){	// regular case
					*(nodal_coeff_pos++) -= interp_quad_2(nodal_pos[-stride + j], nodal_pos[j], nodal_pos[stride + j]);
					*(coeff_coeff_pos++) -= interp_cubic(interp_quad_2(nodal_pos[-stride + j - 1], nodal_pos[j - 1], nodal_pos[stride + j - 1]),
															interp_quad_2(nodal_pos[-stride + j], nodal_pos[j], nodal_pos[stride + j]),
																interp_quad_2(nodal_pos[-stride + j + 1], nodal_pos[j + 1], nodal_pos[stride + j + 1]),
																	interp_quad_2(nodal_pos[-stride + j + 2], nodal_pos[j + 2], nodal_pos[stride + j + 2]));
				}
				// when j = n2_coeff - 1 which means the last coeff element in the last coeff row
				*(nodal_coeff_pos++) -= interp_quad_2(nodal_pos[-stride + n2_coeff - 1], nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1]);
				*(coeff_coeff_pos++) -= interp_quad_2(interp_quad_2(nodal_pos[-stride + n2_coeff - 2], nodal_pos[n2_coeff - 2], nodal_pos[stride + n2_coeff - 2]),
														interp_quad_2(nodal_pos[-stride + n2_coeff - 1], nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1]),
															interp_quad_2(nodal_pos[-stride + n2_coeff], nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff]));
				// compute the last (or second last if n2 is even) nodal column
				*(nodal_coeff_pos++) -= interp_quad_2(nodal_pos[-stride + n2_coeff], nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff]);
				if (even_n2){
					// compute the last nodal column
					*(nodal_coeff_pos++) -= interp_quad_2(nodal_pos[-stride + n2_coeff + 1], nodal_pos[n2_coeff + 1], nodal_pos[stride + n2_coeff + 1]);
				}
			}
			else{	// 0 < i < n1_coeff means regular case
				// when j = 0 means the first coeff element in regular coeff row
				*(nodal_coeff_pos++) -= interp_cubic(nodal_pos[-stride], nodal_pos[0], nodal_pos[stride], nodal_pos[2*stride]);
				*(coeff_coeff_pos++) -= interp_quad_1(interp_cubic(nodal_pos[-stride], nodal_pos[0], nodal_pos[stride], nodal_pos[2*stride]),
														interp_cubic(nodal_pos[-stride + 1], nodal_pos[1], nodal_pos[stride + 1], nodal_pos[2*stride + 1]),
															interp_cubic(nodal_pos[-stride + 2], nodal_pos[2], nodal_pos[stride + 2], nodal_pos[2*stride + 2]));
				for(int j=1; j+1<n2_coeff; j++){
					// coefficients in nodal columns
                	*(nodal_coeff_pos++) -= interp_cubic(nodal_pos[-stride + j], nodal_pos[j], nodal_pos[stride + j], nodal_pos[2*stride + j]);
                	// coefficients in centers
                	*(coeff_coeff_pos++) -= interp_cubic(interp_cubic(nodal_pos[-stride + j - 1], nodal_pos[j - 1], nodal_pos[stride + j - 1], nodal_pos[2*stride + j - 1]),
															interp_cubic(nodal_pos[-stride + j], nodal_pos[j], nodal_pos[stride + j], nodal_pos[2*stride + j]),
																interp_cubic(nodal_pos[-stride + j + 1], nodal_pos[j + 1], nodal_pos[stride + j + 1], nodal_pos[2*stride + j + 1]),
																	interp_cubic(nodal_pos[-stride + j + 2], nodal_pos[j + 2], nodal_pos[stride + j + 2], nodal_pos[2*stride + j + 2]));
				}
				// when j = n2_coeff - 1 which means the last coeff element in regular coeff row
				*(nodal_coeff_pos++) -= interp_cubic(nodal_pos[-stride + n2_coeff - 1], nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1], nodal_pos[2*stride + n2_coeff -1]);
				*(coeff_coeff_pos++) -= interp_quad_2(interp_cubic(nodal_pos[-stride + n2_coeff - 2], nodal_pos[n2_coeff - 2], nodal_pos[stride + n2_coeff - 2], nodal_pos[2*stride + n2_coeff - 2]),
														interp_cubic(nodal_pos[-stride + n2_coeff - 1], nodal_pos[n2_coeff - 1], nodal_pos[stride + n2_coeff - 1], nodal_pos[2*stride + n2_coeff - 1]),
															interp_cubic(nodal_pos[-stride + n2_coeff], nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff], nodal_pos[2*stride + n2_coeff]));
				// compute the last (or second last if n2 is even) nodal column
            	*(nodal_coeff_pos++) -= interp_cubic(nodal_pos[-stride + n2_coeff], nodal_pos[n2_coeff], nodal_pos[stride + n2_coeff], nodal_pos[2*stride + n2_coeff]);
				if (even_n2){
					// compute the last nodal column
					*(nodal_coeff_pos++) -= interp_cubic(nodal_pos[-stride + n2_coeff + 1], nodal_pos[n2_coeff + 1], nodal_pos[stride + n2_coeff + 1], nodal_pos[2*stride + n2_coeff + 1]);
				}
			}
		}
	}
	void compute_interpolant_difference_2D(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		// compute horizontal difference
		const T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n2_nodal;
		for(int i=0; i<n1_nodal; i++){
			compute_interpolant_difference_1D(n2_coeff, nodal_pos, coeff_pos);
			nodal_pos += stride, coeff_pos += stride;
		}
		// compute vertical difference
		compute_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
	}
	void compute_interpolant_difference_2D_cubic(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		// compute horizontal difference
		const T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n2_nodal;
		for(int i=0; i<n1_nodal; i++)
		{
			compute_interpolant_difference_1D_cubic(n2_coeff, nodal_pos, coeff_pos);
			nodal_pos += stride, coeff_pos += stride;
		}
		compute_interpolant_difference_2D_vertical_cubic(data_pos, n1, n2, stride);
	}	
	// decompose n1 x n2 data into coarse level (n1/2 x n2/2)
	void decompose_level_2D(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		// cerr << "decompose, h = " << h << endl; 
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n1_coeff = n1 - n1_nodal;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n2_coeff = n2 - n2_nodal;
		data_reorder_2D(data_pos, data_buffer, n1, n2, stride);
		compute_interpolant_difference_2D(data_pos, n1, n2, stride);
        vector<T> w1(n1_nodal);
        vector<T> b1(n1_nodal);
        vector<T> w2(n2_nodal);
        vector<T> b2(n2_nodal);
        precompute_w_and_b(w1.data(), b1.data(), n1_nodal);
        precompute_w_and_b(w2.data(), b2.data(), n2_nodal);
        compute_correction_2D(data_pos, data_buffer, load_v_buffer, n1, n2, n1_nodal, h, stride, w1.data(), b1.data(), w2.data(), b2.data(), default_batch_size);
        apply_correction_batched(data_pos, data_buffer, n1_nodal, stride, n2_nodal, true);
	}
    void decompose_level_2D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
        data_reorder_2D(data_pos, data_buffer, n1, n2, stride);
        compute_interpolant_difference_2D(data_pos, n1, n2, stride);
    }
	void decompose_level_2D_cubic(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		data_reorder_2D(data_pos, data_buffer, n1, n2, stride);
		compute_interpolant_difference_2D_cubic(data_pos, n1, n2, stride);
	}
	/*
		2D computation + vertical computation for coefficient plane 
	*/
	void compute_interpolant_difference_3D(T * data_pos, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		size_t n3_nodal = (n3 >> 1) + 1;
		size_t n3_coeff = n3 - n3_nodal;
		bool even_n2 = (!(n2 & 1));
		bool even_n3 = (!(n3 & 1));
		T * cur_data_pos = data_pos;
 		for(int i=0; i<n1_nodal; i++){
 			compute_interpolant_difference_2D(cur_data_pos, n2, n3, dim1_stride);
 			cur_data_pos += dim0_stride;
 		}
 		// compute vertically
 		const T * nodal_pos = data_pos;
 		T * coeff_pos = data_pos + n1_nodal * dim0_stride;
		for(int i=0; i<n1_coeff; i++){
			// iterate throught coefficient planes along n1
			/*
				data in the coefficient plane
				xxxxx		xxx						xx
				xxxxx		xxx	coeff_nodal_nonal	xx 	coeff_nodal_coeff
				xxxxx	=>	xxx						xx
				xxxxx
				xxxxx		xxx	coeff_coeff_nodal	xx 	coeff_coeff_coeff
							xxx						xx
			*/
            const T * nodal_nodal_nodal_pos = nodal_pos;
            T * coeff_nodal_nodal_pos = coeff_pos;
            T * coeff_nodal_coeff_pos = coeff_pos + n3_nodal;
            T * coeff_coeff_nodal_pos = coeff_pos + n2_nodal * dim1_stride;
            T * coeff_coeff_coeff_pos = coeff_coeff_nodal_pos + n3_nodal;
            // TODO: optimize average computation
            for(int j=0; j<n2_coeff; j++){
            	for(int k=0; k<n3_coeff; k++){
	                // coeff_nodal_nonal
					//if(i == 50 && j == 50 && k > 0 && k < n3_coeff - 1) std::cout << coeff_nodal_nodal_pos[0] << "\t" << (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2 << std::endl;
	                coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                // coeff_nodal_coeff
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
	                // coeff_coeff_nodal
	                coeff_coeff_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride]) / 4;
	                // coeff_coeff_coeff
	                coeff_coeff_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1] +
	                								nodal_nodal_nodal_pos[k + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride] + 
	                								nodal_nodal_nodal_pos[k + dim1_stride + 1] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride + 1]) / 8;
                }
	            // compute the last (or second last if n3 is even) coeff_*_nodal column
	            coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
                coeff_coeff_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff] +
                								nodal_nodal_nodal_pos[n3_coeff + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + dim1_stride]) / 4;
	            if(even_n3){
	            	// compute the last coeff_*_nodal column if n3 is even
		            coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
	                coeff_coeff_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1] +
	                								nodal_nodal_nodal_pos[n3_coeff + 1 + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1 + dim1_stride]) / 4;
	            }
	            coeff_nodal_nodal_pos += dim1_stride;
	            coeff_nodal_coeff_pos += dim1_stride;
	            coeff_coeff_nodal_pos += dim1_stride;
	            coeff_coeff_coeff_pos += dim1_stride;
	            nodal_nodal_nodal_pos += dim1_stride;
            }
            // compute the last (or second last if n2 is even) coeff_nodal_* row
            {
            	for(int k=0; k<n3_coeff; k++){
            		coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
            	}
        		coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
            	if(even_n3){
            		coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
            	}
                coeff_nodal_nodal_pos += dim1_stride;
                coeff_nodal_coeff_pos += dim1_stride;
                coeff_coeff_nodal_pos += dim1_stride;
                coeff_coeff_coeff_pos += dim1_stride;
                nodal_nodal_nodal_pos += dim1_stride;
        	}
            if(even_n2){
                // compute the last coeff_nodal_* row if n2 is even
            	for(int k=0; k<n3_coeff; k++){
            		coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
            	}
        		coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
            	if(even_n3){
            		coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
            	}
            }
            nodal_pos += dim0_stride;
            coeff_pos += dim0_stride;
		}
	}
	void compute_interpolant_difference_3D_cubic(T * data_pos, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		size_t n3_nodal = (n3 >> 1) + 1;
		size_t n3_coeff = n3 - n3_nodal;
		bool even_n2 = (!(n2 & 1));
		bool even_n3 = (!(n3 & 1));
		T * cur_data_pos = data_pos;
 		for(int i=0; i<n1_nodal; i++){
 			compute_interpolant_difference_2D_cubic(cur_data_pos, n2, n3, dim1_stride);
 			cur_data_pos += dim0_stride;
 		}
 		// compute vertically
 		const T * nodal_pos = data_pos;
 		T * coeff_pos = data_pos + n1_nodal * dim0_stride;
		for(int i=0; i<n1_coeff; i++){
			// iterate throught coefficient planes along n1
			/*
				data in the coefficient plane
				xxxxx		xxx						xx
				xxxxx		xxx	coeff_nodal_nonal	xx 	coeff_nodal_coeff
				xxxxx	=>	xxx						xx
				xxxxx
				xxxxx		xxx	coeff_coeff_nodal	xx 	coeff_coeff_coeff
							xxx						xx
			*/
            const T * nodal_nodal_nodal_pos = nodal_pos;
            T * coeff_nodal_nodal_pos = coeff_pos;
            T * coeff_nodal_coeff_pos = coeff_pos + n3_nodal;
            T * coeff_coeff_nodal_pos = coeff_pos + n2_nodal * dim1_stride;
            T * coeff_coeff_coeff_pos = coeff_coeff_nodal_pos + n3_nodal;
            // TODO: optimize average computation
			if (i == 0){	// which means the first coeff plane																							2**
				for(int j=0; j<n2_coeff; j++){
					if(j == 0){	// which means the first coeff row in the first coeff plane																	211 213 215 217 218 | 212 214 216 && 221 223 225 227 228 | 222 224 226
						// k = 0 which means the first coeff element in the first coeff row in the first coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 211 = q1(111, 311, 511)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																  interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	// 212 = q1(211, 213, 215) = q1(q1(111, 311, 511), q1(113, 313, 513), q1(115, 315, 515))
						coeff_coeff_nodal_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),
																  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride]));	// 221 = q1(211, 231, 251) = q1(q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551)) 
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]), 																// q1(111, 311, 511)
																				interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]), 						// q1(131, 331, 531)	q1(211, 231, 251)	221
																				interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride])), 				// q1(151, 351, 551)
																  interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]), 														// q1(113, 313, 513)
																				interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1]), 			// q1(133, 333, 533)	q1(213, 233, 253)	223
																				interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 1])), 	// q1(153, 353, 553)
																  interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]), 														// q1(115, 315, 515)
																				interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2]), 			// q1(135, 335, 535)	q1(215, 235, 255)	225
																				interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 2]))); 	// q1(155, 355, 555)
																				// 222 = q1(221, 223, 225) = q1(q1(211, 231, 251), q1(213, 233, 253), q1(215, 235, 255)) #9
																				//= q1{q1[q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551)], q1[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553)], q1[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555)]} #27
						for(int k=1; k+1<n3_coeff; k++){	// regular case in the fisrt coeff row in the first coeff plane
							coeff_nodal_nodal_pos[k] -= interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]); // 213 = q1(113, 313, 513)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]), 
																	 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]), 
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 2 ], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2])); // 214 = c(211, 213, 215, 217) = c(q1(111, 311, 511), q1(113, 313, 513), q1(115, 315, 515), q1(117, 317, 517))
							coeff_coeff_nodal_pos[k] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),
																	  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k])); 
																	  // 223 = q1(213, 233, 253) =  q1(q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),													// q1(111, 311, 511)
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1]),			// q1(131, 331, 531)	q1(211, 231, 251)	221
																				   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k - 1])),	// q1(151, 351, 551)
																	 interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),																// q1(113, 313, 513)
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),						// q1(133, 333, 533)	q1(213, 233, 253)	223
																				   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k])),				// q1(153, 353, 553)
																	 interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),													// q1(115, 315, 515)
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1]),			// q1(135, 335, 535)	q1(215, 235, 255)	225
																				   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 1])),	// q1(155, 355, 555)
																	 interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]),													// q1(117, 317, 517)
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2]),			// q1(137, 337, 537)	q1(217, 237, 257)	227
																				   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 2])));	// q1(157, 357, 557)
																	  // 224 = c(221, 223, 225, 227) =  c(q1(211, 231, 251), q1(213, 233, 253), q1(215, 235, 255), q1(217, 237, 257))
																	  //= c{q1[q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551)], q1[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553)], q1[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555)], q1[q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 557)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the first coeff row in the first coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]); // 215 = q1(115, 315, 515)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff])); // 216 = q2(213, 215, 217) = q2(q1(113, 313, 513), q1(115, 315, 515), q1(117, 317, 517))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1]));
																			 // 225 = q1(215, 235, 255) = q1(q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),													// q1(113, 313, 513)	
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2]),		// q1(133, 333, 533)	q1(213, 233, 253)	223
																						   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 2])),	// q1(153, 353, 553)
																			 interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),													// q1(115, 315, 515)
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),		// q1(135, 335, 535)	q1(215, 235, 255)	225
																						   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1])),	// q1(155, 355, 555)
																			 interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),																// q1(117, 317, 517)
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]), 					// q1(137, 337, 537)	q1(217, 237, 257)	227
																						   interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff]))); 			// q1(157, 357, 557)
																			 // 226 = q2(223, 225, 227) = q2(q1(213, 233, 253), q1(215, 235, 255), q1(217, 237, 257))
																			 // = q2{q1[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553)], q1[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555)], q1[q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 557)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]); // 217 = q1(117, 317, 517)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),
																		 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),
																		 interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff])); 
																		 // 227 = q1(217, 237, 257) = q1[q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 557)]
						if(even_n3){ // compute the last coeff_*_nodal column if n3 is even
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]); // 218 = q1(118, 318, 518)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),
																				 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1]),
																				 interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff + 1])); 
																				 // 228 = q1(218, 238, 258) = q1(q1(118, 318, 518), q1(138, 338, 538), q1(158, 358, 558))
						}
					}
					else if(j == n2_coeff-1){	// which means the last coeff row in the first coeff plane													251 253 255 257 258 | 252 254 256 && 261 263 265 267 268 | 262 264 266
						// k = 0 which means the first coeff element in the last coeff row in the first coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]); // 251 = q1(151, 351, 551) 
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2])); // 252 = q1(251, 253, 255) = q1(q1(151, 351, 551), q1(153, 353, 553), q1(155, 355, 555))
						coeff_coeff_nodal_pos[0] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride])); // 261 = q2(231, 251, 271) = q2(q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),				// q1(131, 331, 531)
																			   interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),														// q1(151, 351, 551)	q2(231, 251, 271)	261
																			   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride])),				// q1(171, 371, 571)
																 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 1]),	// q1(133, 333, 533)
																			   interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),												// q1(153, 353, 553)	q2(233, 253, 273)	263
																			   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1])),	// q1(173, 373, 573)
																 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 2]),	// q1(135, 335, 535)
																			   interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]),												// q1(155, 355, 555)	q2(235, 255, 275)	265
																			   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2]))); 	// q1(175, 375, 575)
																			   // 262 = q1(261, 263, 265) = q1(q2(231, 251, 271), q2(233, 253, 273), q2(235, 255, 275))
													 						   // = q1{q2[q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571)], q2[q1(133, 333, 533), q1(153, 353, 553),q1(173, 373, 573)], q2[q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular case in the last coeff row in the first coeff plane		253
							coeff_nodal_nodal_pos[k] -= interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 253 = q1(153, 353, 553)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	// 254 = c(251, 253, 255, 257) = c(q1(151, 351, 551), q1(153, 353, 553), q1(155, 355, 555), q1(157, 357, 557))
							coeff_coeff_nodal_pos[k] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),
																	  interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k])); 
																	  // 263 = q2(233, 253, 273) = q2(q1(133, 333, 533), q1(153, 353, 553),q1(173, 373, 573))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k - 1]),	// q1(131, 331, 531)
																				   interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),												// q1(151, 351, 551)	q2(231, 251, 271)	261
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1])),	// q1(171, 371, 571)
																	 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),				// q1(133, 333, 533)	
																				   interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]), 															// q1(153, 353, 553)	q2(233, 253, 273)	263
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k])),				// q1(173, 373, 573)
																	 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 1]), 	// q1(135, 335, 535)	
																				   interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]), 												// q1(155, 355, 555)	q2(235, 255, 275)	265
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1])),	// q1(175, 375, 575)
																	 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 2]), 	// q1(137, 337, 537)
																				   interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]), 												// q1(157, 357, 557)	q2(237, 257, 277)	267
																				   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2])));	// q1(177, 377, 577)
														  // 264 = c(261, 263, 265, 267) = c(q2(231, 251, 271), q2(233, 253, 273), q2(235, 255, 275), q2(237, 257, 277))
														  // = c{q2[q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571)], q2[q1(133, 333, 533), q1(153, 353, 553),q1(173, 373, 573)], q2[q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)], q2[q1(137, 337, 537), q1(157, 357, 557), q1(177, 377, 577)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the last coeff row in the first coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 255 = q1(155, 355, 555)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	
																			 // 256 = q2(253, 255, 257) = q2(q1(153, 353, 553), q1(155, 355, 555), q1(157, 357, 557))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]));	
																			 // 265 = q2(235, 255, 275) = q2(q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 2]),	// q1(133, 333, 533)
																						   interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),												// q1(153, 353, 553)	q2(233, 253, 273)	263
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2])),	// q1(173, 373, 573)
																			 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),	// q1(135, 335, 535)
																						   interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),												// q1(155, 355, 555)	q2(235, 255, 275)	265
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1])),	// q1(175, 375, 575)
																			 interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),				// q1(137, 337, 537)
																						   interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),															// q1(157, 357, 557)	q2(237, 257, 277)	267
																						   interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]))); 				// q1(177, 377, 577)
																 // 266 = q2(263, 265, 267) = q2(q2(233, 253, 273), q2(235, 255, 275), q2(237, 257, 277))
																 // = q2{q2[q1(133, 333, 533), q1(153, 353, 553), q1(173, 373, 573)], q2[q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)], q2[q1(137, 337, 537), q1(157, 357, 557), q1(177, 377, 577)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 257 = q1(157, 357, 557)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),
																		 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),
																		 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]));	
																		 // 267 = q2(237, 257, 277) = q2(q1(137, 337, 537), q1(157, 357, 557), q1(177, 377, 577))
						if(even_n3){	// compute the last coeff_*_nodal column if n3 is even
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]); // 258 = q1(158, 358, 558)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff + 1]),
																				 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),
																				 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1])); 
																				 // 268 = q2(238, 258, 278) = q2(q1(138, 338, 538), q1(158, 358, 558), q1(178, 378, 578))
						}
					}
					else{	// regular row in the first coeff plane																							231 233 235 237 238 | 232 234 236 && 241 243 245 247 248 | 242 244 246
						// k = 0 which means the first coeff element in the middle coeff row in the first coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 231 = q1(131, 331, 531)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	// 232 = q1(231, 233, 235) = q1(q1(131, 331, 531), q1(133, 333, 533), q1(135, 335, 535))
						coeff_coeff_nodal_pos[0] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),
																 interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride]));	
																 // 241 = c(211, 231, 251, 271) = c(q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),					// q1(111, 311, 511)
																			  interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),															// q1(131, 331, 531)	c(211, 231, 251, 271)	241
																			  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),						// q1(151, 351, 551)
																			  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride])),				// q1(171, 371, 571)	****
																 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 1]),			// q1(113, 313, 513)
																			  interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),													// q1(133, 333, 533)	c(213, 233, 253, 273)	243
																			  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1]),			// q1(153, 353, 553)
																			  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 1])),	// q1(173, 373, 573)	****
																 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 2]),			// q1(115, 315, 515)
																			  interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]),													// q1(135, 335, 535)	c(215, 235, 255, 275)	245
																			  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2]),			// q1(155, 355, 555)
																			  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 2])));	// q1(175, 375, 575)	****
														// 242 = q1(241, 243, 245) = q1(c(211, 231, 251, 271), c(213, 233, 253, 273), c(215, 235, 255, 275))
														// = q1{c[q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571)], c[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553), q1(173, 373, 573)], c[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular case in the middle coeff row in the first coeff plane
							coeff_nodal_nodal_pos[k] -= interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 233 = q1(133, 333, 533)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	// 234 = c(231, 233, 235, 237) = c(q1(131, 331, 531), q1(133, 333, 533), q1(135, 335, 535), q1(137, 337, 537))
							coeff_coeff_nodal_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),
																	 interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k]));	
																	 // 243 = c(213, 233, 253, 273) = c(q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553), q1(173, 373, 573))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k - 1]),			// q1(111, 311, 511)
																				  interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),													// q1(131, 331, 531)	c(211, 231, 251, 271)	241
																				  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1]),			// q1(151, 351, 551)
																				  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k - 1])),	// q1(171, 371, 571)	****
																	 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),						// q1(113, 313, 513)
																				  interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),																// q1(133, 333, 533)	c(213, 233, 253, 273)	243
																				  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),						// q1(153, 353, 553)
																				  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k])),				// q1(173, 373, 573)	****
																	 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 1]),			// q1(115, 315, 515)
																				  interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),													// q1(135, 335, 535)	c(215, 235, 255, 275)	245
																				  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1]),			// q1(155, 355, 555)
																				  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 1])),	// q1(175, 375, 575)	****
																	 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 2]),			// q1(117, 317, 517)
																				  interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]),													// q1(137, 337, 537)	c(217, 237, 257, 277)
																				  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2]),			// q1(157, 357, 557)
																				  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 2])));	// q1(177, 377, 577)	****
														// 244 = c(241, 243, 245, 247) = c(c(211, 231, 251, 271), c(213, 233, 253, 273), c(215, 235, 255, 275), c(217, 237, 257, 277))
														// = c{c[q1(111, 311, 511), q1(131, 331, 531), q1(151, 351, 551), q1(171, 371, 571)], c[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553), q1(173, 373, 573)], c[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)], c[q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 557), q1(177, 377, 577)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the middle coeff row in the first coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 235 = q1(135, 335, 735)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	// 236 = q2(233, 235, 237) = q2(q1(133, 333, 533), q1(135, 335, 535), q1(137, 337, 537))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),
																			interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),
																			interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1]));	
																			// 245 = c(215, 235, 255, 275) = c(q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 2]),		// q1(113, 313, 513)
																						  interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),													// q1(133, 333, 533)	c(213, 233, 253, 273)	243
																						  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2]),			// q1(153, 353, 553)	
																						  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 2])),	// q1(173, 373, 573)	****
																			 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),		// q1(115, 315, 515)
																						  interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),													// q1(135, 335, 535)	c(215, 235, 255, 275)	245
																						  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),			// q1(155, 355, 555)
																						  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1])),	// q1(175, 375, 575)	****
																			 interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),					// q1(117, 317, 517)
																						  interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),																// q1(137, 337, 537)	c(217, 237, 257, 277)	247
																						  interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),						// q1(157, 357, 757)
																						  interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff])));				// q1(177, 377, 577)	****
																			// 246 = q2(243, 245, 247) = q2(c(213, 233, 253, 273), c(215, 235, 255, 275), c(217, 237, 257, 277))
																			// = q2{c[q1(113, 313, 513), q1(133, 333, 533), q1(153, 353, 553), q1(173, 373, 573)], c[q1(115, 315, 515), q1(135, 335, 535), q1(155, 355, 555), q1(175, 375, 575)], c[q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 757), q1(177, 377, 577)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 237 = q1(137, 337, 537)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),
																		interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),
																		interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),
																		interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff]));	
																		// 247 = c(217, 237, 257, 277) = c(q1(117, 317, 517), q1(137, 337, 537), q1(157, 357, 557), q1(177, 377, 577))
						if(even_n3){	// compute the last coeff_*_nodal column if n3 is even
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 238 = q1(138, 338, 538)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff + 1]),
																				interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),
																				interp_quad_1(nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1]),
																				interp_quad_1(nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff + 1]));	
																				// 248 = c(218, 238, 258, 278) = c(q1(118, 318, 518), q1(138, 338, 538), q1(158, 358, 558), q1(178, 378, 578))
						}
					}
					coeff_nodal_nodal_pos += dim1_stride;
	            	coeff_nodal_coeff_pos += dim1_stride;
	            	coeff_coeff_nodal_pos += dim1_stride;
	            	coeff_coeff_coeff_pos += dim1_stride;
	            	nodal_nodal_nodal_pos += dim1_stride;
				}
				// compute the last (or second last if n2 is even) coeff_nodal_* row in the first plane														271 273 275 277 278 | 272 274 276 
				if(1){
					// k = 0 which means the first coeff element in the last (or second last if n2 is even) coeff_nodal_* row in the first plane
					coeff_nodal_nodal_pos[0] -= interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 271 = q1(171, 371 571)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
															  interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
															  interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	// 272 = q1(271, 273, 275) = q1(q1(171, 371, 571), q1(173, 373, 573), q1(175, 375, 575))
					for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the last (or second last if n2 is even) coeff_nodal_* row in the first plane
						coeff_nodal_nodal_pos[k] -= interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 273 = q1(173, 373, 573)
						coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																 interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	// 274 = c(271, 273, 275, 277) = c(q1(171, 371, 571), q1(173, 373, 573), q1(175, 375, 575), q1(177, 377, 577))
					}
					// k = n3_coeff - 1 which means the last coeff elements in the last (or second last if n2 is even) coeff_nodal_* row in the first plane
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 275 = q1(175, 375, 575)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																		interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																		interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	// 276 = q2(273, 275, 277) = q2(q1(173, 373, 573), q1(175, 375, 575), q1(177, 377, 577))
					// compute the last (or second last if n3 is even) coeff_*_nodal column
					coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 277 = q1(177, 377, 577)
					if(even_n3){
						coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 278 = q1(178, 378, 578)
					}
					coeff_nodal_nodal_pos += dim1_stride;
	                coeff_nodal_coeff_pos += dim1_stride;
    	            coeff_coeff_nodal_pos += dim1_stride;
        	        coeff_coeff_coeff_pos += dim1_stride;
            	    nodal_nodal_nodal_pos += dim1_stride;
				}
				if(even_n2){	// compute the last coeff_nodal_* row if n2 is even																			281 283 285 287 288 | 282 284 286
					// k = 0 which means the first element in the last coeff_nodal_* row in the first plane
					coeff_nodal_nodal_pos[0] -= interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 281 = q1(181, 381, 581)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_1(nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
															  interp_quad_1(nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
															  interp_quad_1(nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	// 282 = q1(281, 283, 285) = q1(q1(181, 381, 581), q1(183, 383, 583), q1(185, 385, 585))
            		for(int k=1; k+1<n3_coeff; k++){	// regular middle elements in the last coeff_nodal_* row in the first plane
            			coeff_nodal_nodal_pos[k] -= interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 283 = q1(183, 383, 583)
	                	coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_1(nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																 interp_quad_1(nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																 interp_quad_1(nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	// 284 = c(281, 283, 285, 287) = c(q1(181, 381, 581), q1(183, 383, 583), q1(185, 385, 585), q1(187, 387, 587))
            		}
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 285 = q1(185, 385, 585)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																		 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																		 interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	// 286 = q2(283, 285, 287) = q2(q1(183, 383, 583), q1(185, 385, 585), q1(187, 387, 587))
        			coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 287 = q1(187, 387, 587)
            		if(even_n3){
            			coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_1(nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 288 = q1(188, 388, 588)
            		}
            	}
			}
			else if(i == n1_coeff - 1){	// which means the last coeff plane																					6**
				for(int j=0; j<n2_coeff; j++){
					if(j == 0){	// which means the first coeff row in the last coeff plane																	611 613 615 617 618 | 612 614 616 && 621 623 625 627 628 | 622 624 626
						// k = 0 which means the first coeff element in the first coeff row in the last coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]);	// 611 = q2(311, 511, 711)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]));	// 612 = q1(611, 613, 615) = q1(q2(311, 511, 711), q2(313, 513, 713), q2(315, 515, 715))
						coeff_coeff_nodal_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride]));	
																  // 621 = q1(611, 631, 651) = q1(q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),																// q2(311, 511, 711)
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride]),							// q2(331, 531, 731)	q1(611, 631, 651)	621
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride])),					// q2(351, 551, 751)
																  interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),														// q2(313, 513, 713)
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1]),				// q2(333, 533, 733)	q1(613, 633, 653)	623
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1])),		// q2(353, 553, 753)
																  interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]),														// q2(315, 515, 715)
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2]),				// q2(335, 535, 735)	q1(615, 635, 655)	625
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2])));		// q2(355, 555, 755)
														// 622 = q1(621, 623, 625) = q1(q1(611, 631, 651), q1(613, 633, 653), q1(615, 635, 655))
														// = q1{q1[q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751)], q1[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753)], q1[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the first coeff row in the last coeff plane
							coeff_nodal_nodal_pos[k] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]);	// 613 = q2(313, 513, 713)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]));	
																	 // 614 = c(611, 613, 615, 617) = c(q2(311, 511, 711), q2(313, 513, 713), q2(315, 515, 715), q2(317, 517, 717))
							coeff_coeff_nodal_pos[k] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k]),
																	  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k]));	
																	  // 623 = q1(613, 633, 653) = q1(q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride]),															// q2(311, 511, 711)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1]),			// q2(331, 531, 731)	q1(611, 631, 651)	621
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1])),	// q2(351, 551, 751)
																	 interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),																// q2(313, 513, 713)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k]),						// q2(333, 533, 733)	q1(613, 633, 653)	623
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k])),				// q2(353, 553, 753)
																	 interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),													// q2(315, 515, 715)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1]),			// q2(335, 535, 735)	q1(615, 635, 655)	625
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1])),	// q2(355, 555, 755)
																	 interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]),													// q2(317, 517, 717)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2]),			// q2(337, 537, 737)	q1(617, 637, 657)	627
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2])));	// q2(357, 557, 757)
															// 624 = c(621, 623, 625, 627) = c(q1(611, 631, 651), q1(613, 633, 653), q1(615, 635, 655), q1(617, 637, 657))
															// = c{q1[q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751)], q1[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753)], q1[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755)], q1[q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the first row in the last coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 615 = q2(315, 515, 715)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]));	
																			 // 616 = q2(613, 615, 617) = q2(q2(313, 513, 713), q2(315, 515, 715), q2(317, 517, 717))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1]));	
																			 // 625 = q1(615, 635, 655) = q1(q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),													// q2(313, 513, 713)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2]),			// q2(333, 533, 733)	q1(613, 633, 653)	623
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2])),	// q2(353, 553, 753)
																			 interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),													// q2(315, 515, 715)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1]),			// q2(335, 535, 735)	q1(615, 635, 655)	625
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1])),	// q2(355, 555, 755)
																			 interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),																// q2(317, 517, 717)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff]),						// q2(337, 537, 737)	q1(617, 637, 657)	627
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff])));				// q2(357, 557, 757)
																	// 626 = q2(623, 625, 627) = q2(q1(613, 633, 653), q1(615, 635, 655), q1(617, 637, 657))
																	// = q2{q1[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753)], q1[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755)], q1[q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]);	// 617 = q2(317, 517, 717)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff]));	
																		 // 627 = q1(617, 637, 657) = q1(q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757))
						if(even_n3){	// compute the last coeff_*_nodal column if n3 is even
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]);	// 618 = q2(318, 518, 718)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]),
																				 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1]),
																				 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff + 1]));	
																				 // 628 = q1(618, 638, 658) = q1(q2(318, 518, 718), q2(338, 538, 738), q2(358, 558, 758))
						}
					}
					else if (j == n2_coeff-1){	// which means the last coeff row in the last coeff plane													651 653 655 657 658 | 652 654 656 && 661 663 665 667 668 | 662 664 666
						// k = 0 which means the first coeff element in the last coeff row in the last coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]);	// 651 = q2(351, 551, 751)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]));
																  // 652 = q1(651, 653, 655) = q1(q2(351, 551, 751), q2(353, 553, 753), q2(355, 555, 755))
						coeff_coeff_nodal_pos[0] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride]));
																  // 661 = q2(631, 651, 671) = q2(q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride]),					// q2(331, 531, 731)
																  				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),															// q2(351, 551, 751)	q2(631, 651, 671)	661
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride])),					// q2(371, 571, 771)
																  interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1]),		// q2(333, 533, 733)
																  				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),													// q2(353, 553, 753)	q2(633, 653, 673)	663
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1])),		// q2(373, 573, 773)
																  interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2]),		// q2(335, 535, 735)
																  				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]),													// q2(355, 555, 755)	q2(635, 655, 675)	665
																				interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2])));		// q2(375, 575, 775)
																	// 662 = q1(661, 663, 665) = q1(q2(631, 651, 671), q2(633, 653, 673), q2(635, 655, 675))
																	// = q1{q2[q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771)], q2[q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)], q2[q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the last ceoff row in the last coeff plane
							coeff_nodal_nodal_pos[k] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]);	// 653 = q2(353, 553, 753)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]));	
																	 // 654 = c(651, 653, 655, 657) = c(q2(351, 551, 751), q2(353, 553, 753), q2(355, 555, 755), q2(357, 557, 757))
							coeff_coeff_nodal_pos[k] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k]),
																	  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k]));
																	  	// 663 = q2(633, 653, 673) = q2(q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1]),		// q2(331, 531, 731)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),												// q2(351, 551, 751)	q2(631, 651, 671)	661
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1])),		// q2(371, 571, 771)
																	 interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k]),					// q2(333, 533, 733)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),															// q2(353, 553, 753)	q2(633, 653, 673)	663
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k])),					// q2(373, 573, 773)
																	 interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1]),		// q2(335, 535, 735)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),												// q2(355, 555, 755)	q2(635, 655, 675)	665
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1])),		// q2(375, 575, 775)
																	 interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2]),		// q2(337, 537, 737)
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]),												// q2(357, 557, 757)	q2(637, 657, 677)	667
																				   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2])));	// q2(377, 577, 777)
																	// 664 = c(661, 663, 665, 667) = c(q2(631, 651, 671), q2(633, 653, 673), q2(635, 655, 675), q2(637, 657, 677))
																	// = c{q2[q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771)], q2[q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)], q2[q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)], q2[q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the last coeff row in the last coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 655 = q2(355, 555, 755)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]));	
																			 // 656 = q2(653, 655, 657) = q2(q2(353, 553, 753), q2(355, 555, 755), q2(357, 557, 757))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1]));	
																			 // 665 = q2(635, 655, 675) = q2(q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2]),		// q2(333, 533, 733)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),													// q2(353, 553, 753)	q2(633, 653, 673)	663
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2])),		// q2(373, 573, 773)
																			 interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1]),		// q2(335, 535, 735)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),													// q2(355, 555, 755)	q2(635, 655, 675)	665
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1])),		// q2(375, 575, 775)
																			 interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff]),					// q2(337, 537, 737)
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),																// q2(357, 557, 757)	q2(637, 657, 677)	667
																						   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff])));					// q2(377, 577, 777)
																	// 666 = q2(663, 665, 667) = q2(q2(633, 653, 673), q2(635, 655, 675), q2(637, 657, 677))
																	// = q2{q2[q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)], q2[q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)], q2[q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]);	// 657 = q2(357, 557, 757)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff]));	
																		 // 667 = q2(637, 657, 677) = q2(q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777))
						if(even_n3){
							coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]);	// 658 = q2(358, 558, 758)
							coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1]));	
																			 // 668 = q2(638, 658, 678) = q2(q2(338, 538, 738), q2(358, 558, 758), q2(378, 578, 778))
						}
					}
					else{	// regular row in the first coeff plane																							631 633 635 637 638 | 632 634 636 && 641 643 645 647 648 | 642 644 646
						// k = 0 which means the first coeff element in the regular middle rows in the last coeff plane
						coeff_nodal_nodal_pos[0] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]);	// 631 = q2(331, 531, 731)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),
																  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]));	
																 // 632 = q1(631, 633, 635) = q1(q2(331, 531, 731), q2(333, 533, 733), q2(335, 535, 735))
						coeff_coeff_nodal_pos[0] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride]));	
																 // 641 = c(611, 631, 651, 671) = c(q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride]),							// q2(311, 511, 711)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),																// q2(331, 531, 731)	c(611, 631, 651, 671)	641
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride]),							// q2(351, 551, 751)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride])),					// q2(371, 571, 771)	****
																  interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1]),				// q2(313, 513, 713)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),														// q2(333, 533, 733)	c(613, 633, 653, 673)	643
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1]),				// q2(353, 553, 753)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1])),		// q2(373, 573, 773)	****
																  interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2]),				// q2(315, 515, 715)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]),														// q2(335, 535, 735)	c(615, 635, 655, 675)	645
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2]),				// q2(355, 555, 755)
																			   interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2])));		// q2(375, 575, 775)	****
														// 642 = q1(641, 643, 645) = q1(c(611, 631, 651, 671), c(613, 633, 653, 673), c(615, 635, 655, 675))
														// = q1{c[q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771)], c[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)], c[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the regular middle rows in the last coeff plane
							coeff_nodal_nodal_pos[k] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]);	// 633 = q2(333, 533, 733)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]));	
																	 // 634 = c(631, 633, 635, 637) = c(q2(331, 531, 731), q2(333, 533, 733), q2(335, 535, 735), q2(337, 537, 737))
							coeff_coeff_nodal_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k]),
																	 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k]));	
																	 // 643 = c(613, 633, 653, 673) = c(q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1]),			// q2(311, 511, 711)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),														// q2(331, 531, 731)	c(611, 631, 651, 671)	641
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1]),			// q2(351, 551, 751)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1])),	// q2(371, 571, 771)	****
																	 interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k]),						// q2(313, 513, 713)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),																	// q2(333, 533, 733)	c(613, 633, 653, 673)	643
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k]),						// q2(353, 553, 753)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k])),				// q2(373, 573, 773)	****
																	 interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1]),			// q2(315, 515, 715)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),														// q2(335, 535, 735)	c(615, 635, 655, 675)	645
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1]),			// q2(355, 555, 755)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1])),	// q2(375, 575, 775)	****
																	 interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2]),			// q2(317, 517, 717)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride]),																		// q2(337, 537, 737)	c(617, 637, 657, 677)	647
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2]),			// q2(357, 557, 757)
																				  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2])));	// q2(377, 577, 777)	****
														// 644 = c(641, 643, 645, 647) = c(c(611, 631, 651, 671), c(613, 633, 653, 673), c(615, 635, 655, 675), c(617, 637, 657, 677))
														// = c{c[q2(311, 511, 711), q2(331, 531, 731), q2(351, 551, 751), q2(371, 571, 771)], c[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)],
														//		 c[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)], c[q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the regular middle coeff row in the last coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 635 = q2(335, 535, 735)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]));	
																			 // 636 = q2(633, 635, 637) = q2(q2(333, 533, 733), q2(335, 535, 735), q2(337, 537, 737))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1]));	
																			// 645 = c(615, 635, 655, 675) = c(q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2]),			// q2(313, 513, 713)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),													// q2(333, 533, 733)	c(613, 633, 653, 673)	643
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2]),			// q2(353, 553, 753)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2])),	// q2(373, 573, 773)	****
																			 interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1]),			// q2(315, 515, 715)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),													// q2(335, 535, 735)	c(615, 635, 655, 675)	645
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1]),			// q2(355, 555, 755)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1])),	// q2(375, 575, 775)	****
																			 interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff]),						// q2(317, 517, 717)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),																// q2(337, 537, 737)	c(617, 637, 657, 677)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff]),						// q2(357, 557, 757)
																						  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff])));				// q2(377, 577, 777)	****
														// 646 = q2(643, 645, 647) = q2(c(613, 633, 653, 673), c(615, 635, 655, 675), c(617, 637, 657, 677))
														// = q2{c[q2(313, 513, 713), q2(333, 533, 733), q2(353, 553, 753), q2(373, 573, 773)], c[q2(315, 515, 715), q2(335, 535, 735), q2(355, 555, 755), q2(375, 575, 775)], c[q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]);	// 637 = q2(337, 537, 737)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff]),
																		interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]),
																		interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff]),
																		interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff]));	
																		// 647 = c(617, 637, 657, 677) = c(q2(317, 517, 717), q2(337, 537, 737), q2(357, 557, 757), q2(377, 577, 777))
						if(even_n3){
							coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]);	// 638 = q2(338, 538, 738)
							coeff_coeff_nodal_pos[n3_coeff] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1]),
																			interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff + 1]));	
																			// 648 = c(618, 638, 658, 678) = c(q2(318, 518, 718), q2(338, 538, 738), q2(358, 558, 758), q2(378, 578, 778))
						}
					}
					coeff_nodal_nodal_pos += dim1_stride;
	            	coeff_nodal_coeff_pos += dim1_stride;
	            	coeff_coeff_nodal_pos += dim1_stride;
	            	coeff_coeff_coeff_pos += dim1_stride;
	            	nodal_nodal_nodal_pos += dim1_stride;
				}
				// compute the last (or second last if n2 is even) coeff_nodal_* row in the last plane														671 673 675 677 678 | 672 674 676
				if(1){
					// k = 0 which means the first coeff element in the last (or second last if n2 is even) coeff_nodal_* row in the last plane
					coeff_nodal_nodal_pos[0] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]);	// 671 = q2(371, 571, 771)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
															  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),
															  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]));	
															  // 672 = q1(671, 673, 675) = q1(q2(371, 571, 771), q2(373, 573, 773), q2(375, 575, 775))
					for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the last (or second last if n2 is even) coeff_nodal_* row in the last plane
						coeff_nodal_nodal_pos[k] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]);	// 673 = q2(373, 573, 773)
						coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2])); 
																 // 674 = c(671, 673, 675, 677) = c(q2(371, 571, 771), q2(373, 573, 773), q2(375, 575, 775), q2(377, 577, 777))
					}
					// k = n3_coeff - 1 which means the last coeff elements in the last (or second last if n2 is even) coeff_nodal_* row in the last plane
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 675 = q2(375, 575, 775)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																		 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]));	
																		 // 676 = q2(673, 675, 677) = q2(q2(373, 573, 773), q2(375, 575, 775), q2(377, 577, 777))
					// compute the last (or second last if n3 is even) coeff_*_nodal column
					coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]);	// 677 = q2(377, 577, 777)
					if(even_n3){
						coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]);	// 678 = q2(378, 578, 778)
					}
					coeff_nodal_nodal_pos += dim1_stride;
	                coeff_nodal_coeff_pos += dim1_stride;
    	            coeff_coeff_nodal_pos += dim1_stride;
        	        coeff_coeff_coeff_pos += dim1_stride;
            	    nodal_nodal_nodal_pos += dim1_stride;
				}
				if(even_n2){	// compute the last coeff_nodal_* row if n2 is even																			681 683 685 687 688 | 682 684 686
					// k = 0 which means the first element in the last coeff_nodal_* row in the last plane
					coeff_nodal_nodal_pos[0] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]);	// 681 = q2(381, 581, 781)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride]),
															  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1]),
															  interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2]));	
															  // 682 = q1(681, 683, 685) = q1(q2(381, 581, 781), q2(383, 583, 783), q2(385, 585, 785))
					for(int k=1; k+1<n3_coeff; k++){	// regular middle elements in the last coeff_nodal_* row in the last plane
						coeff_nodal_nodal_pos[k] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]);	// 683 = q2(383, 583, 783)
						coeff_nodal_coeff_pos[k] -= interp_cubic(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1]),
																 interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2]));	
																 // 684 = c(681, 683, 685, 687) = c(q2(381, 581, 781), q2(383, 583, 783), q2(385, 585, 785), q2(387, 587, 787))
					}
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]);	// 685 = q2(385, 585, 785)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2]),
																		interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1]),
																		interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]));
																		// 686 = q2(683, 685, 687) = q2(q2(383, 583, 783), q2(385, 585, 785), q2(387, 587, 787))
					coeff_nodal_nodal_pos[n3_coeff] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff]);	// 687 = q2(387, 587, 787)
					if(even_n3){
						coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_quad_2(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]);	// 688 = q2(388, 588, 788)
					}
				}
			}
			else{	// which means the middle coeff plane																									4**
				for(int j=0; j<n2_coeff; j++){
					if(j == 0){	// which means the first coeff row in the middle coeff plane																411 413 415 417 418 | 412 414 416 && 421 423 425 427 428 | 422 424 426
						// k = 0 which means the first coeff element in the first coeff row in the middle coeff plane
						coeff_nodal_nodal_pos[0] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 411 = c(111, 311, 511, 711)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	
																  // 412 = q1(411, 413, 415) = q1(c(111, 311, 511, 711), c(113, 313, 513, 713), c(115, 315, 515, 715))
						coeff_coeff_nodal_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride]));	
																  // 421 = q1(411, 431, 451) = q1(c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),																					// c(111, 311, 511, 711)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),								// c(131, 331, 531, 731)	q1(411, 431, 451)	421
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride])),						// c(151, 351, 551, 751)
																  interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),																		// c(113, 313, 513, 713)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1]),				// c(133, 333, 533, 733)	q1(413, 433, 453)	423
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 1])),		// c(153, 353, 553, 753)
																  interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]),																		// c(115, 315, 515, 715)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2]),				// c(135, 335, 535, 735)	q1(415, 435, 455)	425
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 2])));	// c(155, 355, 555, 755)
														// 422 = q1(421, 423, 425) = q1(q1(411, 431, 451), q1(413, 433, 453), q1(415, 435, 455))
														// = q1{q1[c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751)], q1[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753)], q1[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the first coeff row in the middle coeff plane
							coeff_nodal_nodal_pos[k] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 413 = c(113, 313, 513, 713)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	
																	 // 414 = c(411, 413, 415, 417) = c(c(111, 311, 511, 711), c(113, 313, 513, 713), c(115, 315, 515, 715), c(117, 317, 517, 717))
							coeff_coeff_nodal_pos[k] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),
																	  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k]));	
																	  // 423 = q1(413, 433, 453) = q1(c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),																				// c(111, 311, 511, 711)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1]),				// c(131, 331, 531, 731)	q1(411, 431, 451)	421
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k - 1])),		// c(151, 351, 551, 751)
																	 interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),																						// c(113, 313, 513, 713)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),								// c(133, 333, 533, 733)	q1(413, 433, 453)	423
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k])),						// c(153, 353, 553, 753)
																	 interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),																		// c(115, 315, 515, 715)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1]),				// c(135, 335, 535, 735)	q1(415, 435, 455)	425
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 1])),		// c(155, 355, 555, 755)
																	 interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]),																		// c(117, 317, 517, 717)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2]),				// c(137, 337, 537, 737)	q1(417, 437, 457)	427
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 2])));		// c(157, 357, 557, 757)
																		// 424 = c(421, 423, 425, 427) = c(q1(411, 431, 451), q1(413, 433, 453), q1(415, 435, 455), q1(417, 437, 457))
																		// = c{q1[c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751)], q1[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753)], 
																		//		q1[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755)], q1[c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the first coeff row in the middle plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]); // 415 = c(115, 315, 515, 715)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]));	
																			 // 416 = q2(413, 415, 417) = q2(c(113, 313, 513, 713), c(115, 315, 515, 715), c(117, 317, 517, 717))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1]));	
																			 // 425 = q1(415, 435, 455) = q1(c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),		// c(113, 313, 513, 713)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2]),		// c(133, 333, 533, 733)	q1(413, 433, 453)	423
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 2])),		// c(153, 353, 553, 753)
																			 interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),		// c(115, 315, 515, 715)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),		// c(135, 335, 535, 735)	q1(415, 435, 455)	425
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1])),		// c(155, 355, 555, 755)
																			 interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),		// c(117, 317, 517, 717)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2]),		// c(137, 337, 537, 737)	q1(417, 437, 457)	427
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 2])));		// c(157, 357, 557, 757)
																	// 426 = q2(423, 425, 427) = q2(q1(413, 433, 453), q1(415, 435, 455), q1(417, 437, 457))
																	// = q2{q1[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753)], q1[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755)], q1[c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 417 = c(117, 317, 517, 717)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff]));	
																		 // 427 = q1(417, 437, 457) = q1(c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757))
						if(even_n3){
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]); // 418 = c(118, 318, 518, 718)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff + 1]));	
																		 // 428 = q1(418, 438, 458) = q1(c(118, 318, 518, 718), c(138, 338, 538, 738), c(158, 358, 558, 758))
						}
					}
					else if (j == n2_coeff-1){	// which means the last coeff row in the middle coeff plane													451 453 455 457 458 | 452 454 456 && 461 463 465 467 468 | 462 464 466
						// k = 0 which means the first coeff element in the last coeff row in the middle coeff plane
						coeff_nodal_nodal_pos[0] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 451 = c(151, 351, 551, 751)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	
																  // 452 = q1(451, 453, 455) = q1(c(151, 351, 551, 751), c(153, 353, 553, 753), c(155, 355, 555, 755))
						coeff_coeff_nodal_pos[0] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[dim0_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]));	
																  // 461 = q2(431, 451, 471) = q2(c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),						// c(131, 331, 531, 731)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),																			// c(151, 351, 551, 751)	q2(431, 451, 471)	461
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride])),						// c(171, 371, 571, 771)
																  interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride+ 1], nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 1]),		// c(133, 333, 533, 733)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),																// c(153, 353, 553, 753)	q2(433, 453, 473)	463
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1])),		// c(173, 373, 573, 773)
																  interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 2]),		// c(135, 335, 535, 735)
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]),																// c(155, 355, 555, 755)	q2(435, 455, 475)	465
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2])));	// c(175, 375, 575, 775)
														// 462 = q1(461, 463, 465) = q1(q2(431, 451, 471), q2(433, 453, 473), q2(435, 455, 475))
														// = q1{q2[c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771)], q2[c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)], q2[c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the last coeff row in the middle coeff plane
							coeff_nodal_nodal_pos[k] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 453 = c(153, 353, 553, 753)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	
																	 // 454 = c(451, 453, 455, 457) = c(c(151, 351, 551, 751), c(153, 353, 553, 753), c(155, 355, 555, 755), c(157, 357, 557, 757))
							coeff_coeff_nodal_pos[k] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),
																	  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																	  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]));	
																	  // 463 = q2(433, 453, 573) = q2(c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773))
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k - 1]),		// c(131, 331, 531, 731)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),																// c(151, 351, 551, 751)	q2(431, 451, 471)	461
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1])),		// c(171, 371, 571, 771)
																	 interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),						// c(133, 333, 533, 733)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),																				// c(153, 353, 553, 753)	q2(433, 453, 473)	463
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k])),						// c(173, 373, 573, 773)
																	 interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 1]),		// c(135, 335, 535, 735)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),																// c(155, 355, 555, 755)	q2(435, 455, 475)	465
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1])),		// c(175, 375, 575, 775)
																	 interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 2]),		// c(137, 337, 537, 737)
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]),																// c(157, 357, 557, 757)	q2(437, 457, 477))	467
																				   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2])));		// c(177, 377, 577, 777)
															// 464 = c(461, 463, 465, 467) = c(q2(431, 451, 471), q2(433, 453, 473), q2(435, 455, 475), q2(437, 457, 477))
															// = c{q2[c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771)], q2[c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)], 
															//		q2[c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)], q2[c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777)]}
						}
						// k = n3_coeff - 1 which means the last coeff element in the last coeff row in the middle coeff plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 455 = c(155, 355, 555, 755)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	
																			 // 456 = q2(453, 455, 457) = q2(c(153, 353, 553, 753), c(155, 355, 555, 755), c(157, 357, 557, 757))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]));	
																			 // 465 = q2(435, 455, 475) = q2(c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 2]),		// c(133, 333, 533, 733)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),																// c(153, 353, 553, 753)	q2(433, 453, 473)	463
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2])),		// c(173, 373, 573, 773)
																			 interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),		// c(135, 335, 535, 735)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),																// c(155, 355, 555, 755)	q2(435, 455, 475)	465
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1])),		// c(175, 375, 575, 775)
																			 interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),						// c(137, 337, 537, 737)
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),																	// c(157, 357, 557, 757)	q2(437, 457, 477)	467
																						   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff])));						// c(177, 377, 577, 777)
															// 466 = q2(463, 465, 467) = q2(q2(433, 453, 473), q2(435, 455, 475), q2(437, 457, 477))
															// = q2{q2[c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)], q2[c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)], q2[c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 457 = c(157, 357, 557, 757)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]));
																		 // 467 = q2(437, 457, 477) = q2(c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777))
						if(even_n3){
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 458 = c(158, 358, 558, 758)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff + 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),
																			 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1]));
																			 // 468 = q2(438, 458, 478) = q2(c(138, 338, 538, 738), c(158, 358, 558, 758), c(178, 378, 578, 778))
						}
					}
					else{	// which means the regular middle coeff row in the middle coeff plane															431 433 435 437 438 | 432 434 436 && 441 443 445 447 448 | 442 444 446
						// k = 0 which means the first coeff element in the regular middle coeff row in the middle coeff plane
						coeff_nodal_nodal_pos[0] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 431 = c(131, 331, 531, 731)
						coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
																  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	
																  // 432 = q1(431, 433, 435) = q1(c(131, 331, 531, 731), c(133, 333, 533, 733), c(135, 335, 535, 735))
						coeff_coeff_nodal_pos[0] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride]));	
																 // 441 = c(411, 431, 451, 471) = c(c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771))
						coeff_coeff_coeff_pos[0] -= interp_quad_1(interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride], nodal_nodal_nodal_pos[-dim1_stride], nodal_nodal_nodal_pos[dim0_stride - dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride]),							// c(111, 311, 511, 711)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),																				// c(131, 331, 531, 731)	c(411, 431, 451, 471)	441
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride], nodal_nodal_nodal_pos[dim1_stride], nodal_nodal_nodal_pos[dim0_stride + dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride]),							// c(151, 351, 551, 751)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim1_stride], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride])),					// c(171, 371, 571, 771)	****
																  interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[-dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 1]),			// c(113, 313, 513, 713)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),																	// c(133, 333, 533, 733)	c(413, 433, 453, 473)	443
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 1]),			// c(153, 353, 553, 753)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim1_stride + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 1])),	// c(173, 373, 573, 773)	****
																  interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[-dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + 2]),			// c(115, 315, 515, 715)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]),																	// c(135, 335, 535, 735)	c(415, 435, 455, 475)	445
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + 2]),			// c(155, 355, 555, 755)
																			   interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim1_stride + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + 2])));	// c(175, 375, 575, 775)	****
														// 442 = q1(441, 443, 445) = q1(c(411, 431, 451, 471), c(413, 433, 453, 473), c(415, 435, 455, 475))
														// = q1{c[c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771)], 
														//		c[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)], 
														//		c[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)]}
						for(int k=1; k+1<n3_coeff; k++){	// regular middle coeff elements in the regular middle coeff row in the middle coeff plane
							//if(i == 50 && j == 50) std::cout <<  coeff_nodal_nodal_pos[k] << "\t" << interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]) << std::endl;
							coeff_nodal_nodal_pos[k] -=	interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 433 = c(133, 333, 533, 733)
							coeff_nodal_coeff_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	
																	 // 434 = c(431, 433, 435, 437) = c(c(131, 331, 531, 731), c(133, 333, 533, 733), c(135, 335, 535, 735), c(137, 337, 537, 737))
							coeff_coeff_nodal_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),
																	 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k]));	
																	 // 443 = c(413, 433, 453, 473) = c(c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773))-dim1_stride + k - 1
							coeff_coeff_coeff_pos[k] -= interp_cubic(interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[-dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k - 1]),			// c(111, 311, 511, 711)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),																	// c(131, 331, 531, 731)	c(411, 431, 451, 471)	441
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k - 1]),			// c(151, 351, 551, 751)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim1_stride + k - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k - 1])),	// c(171, 371, 571, 771)	***
																	 interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[-dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k]),							// c(113, 313, 513, 713)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),																					// c(133, 333, 533, 733)	c(413, 433, 453, 473)	443
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k]),							// c(153, 353, 553, 753)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim1_stride + k], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k])),					// c(173, 373, 573, 773)	****
																	 interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[-dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 1]),			// c(115, 315, 515, 715)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),																	// c(135, 335, 535, 735)	c(415, 435, 455, 475)	445
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 1]),			// c(155, 355, 555, 755)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim1_stride + k + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 1])),	// c(175, 375, 575, 775)	****
																	 interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[-dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + k + 2]),			// c(117, 317, 517, 717)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]),																	// c(137, 337, 537, 737)	c(417, 437, 457, 477)	447
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + k + 2]),			// c(157, 357, 557, 757)
																				  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim1_stride + k + 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + k + 2])));	// c(177, 377, 577, 777)	****
														// 444 = c(441, 443, 445, 447) = c(c(411, 431, 451, 471), c(413, 433, 453, 473), c(415, 435, 455, 475), c(417, 437, 457, 477))
														// = c{c[c(111, 311, 511, 711), c(131, 331, 531, 731), c(151, 351, 551, 751), c(171, 371, 571, 771)], c[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)],
														//		c[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)], c[c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777)]}											
						}
						// k = n3_coeff - 1 which means the last coeff element in the regular middle row in the middle plane
						coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 435 = c(135, 335, 535, 735)
						coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																			interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));
																			// 436 = q2(433, 435, 437) = c(c(133, 333, 533, 733), c(135, 335, 535, 735), c(137, 337, 537, 737))
						coeff_coeff_nodal_pos[n3_coeff - 1] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),
																			interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																			interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),
																			interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1]));	
																			// 445 = c(415, 435, 455, 475) = c(c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775))
						coeff_coeff_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 2]),			// c(113, 313, 513, 713)
																			 			  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),																	// c(133, 333, 533, 733)	c(413, 433, 453, 473)	443
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 2]),			// c(153, 353, 553, 753)
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 2])),	// c(173, 373, 573, 773)	****
																			 interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff - 1]),			// c(115, 315, 515, 715)
																			 			  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),																	// c(135, 335, 535, 735)	c(415, 435, 455, 475)	445
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff - 1]),			// c(155, 355, 555, 755)
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff - 1])),	// c(175, 375, 575, 775)	****
																			 interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),							// c(117, 317, 517, 717)
																			 			  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),																					// c(137, 337, 537, 737)	c(417, 437, 457, 477)	447
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),							// c(157, 357, 557, 757)
																						  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff])));					// c(177, 377, 577, 777)	****
																// 446 = q2(443, 445, 447) = q2(c(413, 433, 453, 473), c(415, 435, 455, 475), c(417, 437, 457, 477))
																// = q2{c[c(113, 313, 513, 713), c(133, 333, 533, 733), c(153, 353, 553, 753), c(173, 373, 573, 773)],
																//		c[c(115, 315, 515, 715), c(135, 335, 535, 735), c(155, 355, 555, 755), c(175, 375, 575, 775)], 
																//		c[c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777)]}
						// compute the last (or second last if n3 is even) coeff_*_nodal column
						coeff_nodal_nodal_pos[n3_coeff] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 437 = c(137, 337, 537, 737)
						coeff_coeff_nodal_pos[n3_coeff] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff]),
																		interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]),																					
																		interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff]),							
																		interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff]));
																		// 447 = c(417, 437, 457, 477) = c(c(117, 317, 517, 717), c(137, 337, 537, 737), c(157, 357, 557, 757), c(177, 377, 577, 777))
						if(even_n3){
							coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 438 = c(1378, 338, 538, 738)
							coeff_coeff_nodal_pos[n3_coeff + 1] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[-dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride - dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride - dim1_stride + n3_coeff + 1]),	
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]),																					
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + dim1_stride + n3_coeff + 1]),							
																				interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + 2*dim1_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + 2*dim1_stride + n3_coeff + 1]));
																				// 448 = c(418, 438, 458, 478) = c(c(118, 318, 518, 718), c(138, 338, 538, 738), c(158, 358, 558, 758), c(178, 378, 578, 778))
						}
					}
					coeff_nodal_nodal_pos += dim1_stride;
	            	coeff_nodal_coeff_pos += dim1_stride;
	            	coeff_coeff_nodal_pos += dim1_stride;
	            	coeff_coeff_coeff_pos += dim1_stride;
	            	nodal_nodal_nodal_pos += dim1_stride;
				}
				// compute the last (or second last if n2 is even) coeff_nodal_* row in the middle plane													471 473 475 477 478 | 472 474 476
				if(1){
					// k = 0 which means the first coeff element in the last (or second last if n2 is even) coeff_nodal_* row in the middle plane
					coeff_nodal_nodal_pos[0] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 471 = c(171, 371, 571, 771)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
															  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
															  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));	
															  // 472 = q1(471, 473, 475) = q1(c(171, 371, 571, 771), c(173, 373, 573, 773), c(175, 375, 575, 775))
					for(int k=1; k+1<n3_coeff; k++){
						coeff_nodal_nodal_pos[k] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 473 = c(173, 373, 573, 773)
						coeff_nodal_coeff_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	
																 // 474 = c(471, 473, 475, 477) = c(c(171, 371, 571, 771), c(173, 373, 573, 773), c(175, 375, 575, 775), c(177, 377, 577, 777))
					}
					// k = n3_coeff - 1 which means the last coeff element in the last (or second last if n2 is even) coeff_nodal_* row in the middle plane
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 475 = c(175, 375, 575, 775)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	
																		 // 476 = q2(473, 475, 477) = q2(c(173, 373, 573, 773), c(175, 375, 575, 775), c(177, 377, 577, 777))
					coeff_nodal_nodal_pos[n3_coeff] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 477 = c(177, 377, 577, 777)
					if(even_n3){
						coeff_nodal_coeff_pos[n3_coeff + 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 478 = c(178, 378, 578, 778)
					}
					coeff_nodal_nodal_pos += dim1_stride;
	            	coeff_nodal_coeff_pos += dim1_stride;
	            	coeff_coeff_nodal_pos += dim1_stride;
	            	coeff_coeff_coeff_pos += dim1_stride;
	            	nodal_nodal_nodal_pos += dim1_stride;
				}
				if(even_n2){	// compute the last coeff_nodal_* row if n2 is even																			481 483 485 487 488 | 482 484 486
					// k = 0 which means the first element in the last coeff_nodal_* row in the middle plane
					coeff_nodal_nodal_pos[0] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]);	// 481 = c(181, 381, 581, 781)
					coeff_nodal_coeff_pos[0] -= interp_quad_1(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride], nodal_nodal_nodal_pos[0], nodal_nodal_nodal_pos[dim0_stride], nodal_nodal_nodal_pos[2*dim0_stride]),
															  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 1], nodal_nodal_nodal_pos[1], nodal_nodal_nodal_pos[dim0_stride + 1], nodal_nodal_nodal_pos[2*dim0_stride + 1]),
															  interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + 2], nodal_nodal_nodal_pos[2], nodal_nodal_nodal_pos[dim0_stride + 2], nodal_nodal_nodal_pos[2*dim0_stride + 2]));
															  // 482 = q1(481, 483, 485) = q1(c(181, 381, 581, 781), c(183, 383, 583, 783), c(185, 385, 585, 785))
					for(int k=1; k+1<n3_coeff; k++){	// regular middle elements in the last coeff_nodal_* row in the middle plane
						coeff_nodal_nodal_pos[k] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]);	// 483 = c(183, 383, 583, 783)
						coeff_nodal_coeff_pos[k] -= interp_cubic(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k - 1], nodal_nodal_nodal_pos[k - 1], nodal_nodal_nodal_pos[dim0_stride + k - 1], nodal_nodal_nodal_pos[2*dim0_stride + k - 1]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k], nodal_nodal_nodal_pos[k], nodal_nodal_nodal_pos[dim0_stride + k], nodal_nodal_nodal_pos[2*dim0_stride + k]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 1], nodal_nodal_nodal_pos[k + 1], nodal_nodal_nodal_pos[dim0_stride + k + 1], nodal_nodal_nodal_pos[2*dim0_stride + k + 1]),
																 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + k + 2], nodal_nodal_nodal_pos[k + 2], nodal_nodal_nodal_pos[dim0_stride + k + 2], nodal_nodal_nodal_pos[2*dim0_stride + k + 2]));	
																 // 484 = c(481, 483, 485, 487) = c(c(181, 381, 581, 781), c(183, 383, 583, 783), c(185, 385, 585, 785), c(187, 387, 587, 787))
					}
					// k = n3_coeff - 1 which means the last element in the last coeff_nodal_* row in the middle plane
					coeff_nodal_nodal_pos[n3_coeff - 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]);	// 485 = c(185, 385, 585, 785)
					coeff_nodal_coeff_pos[n3_coeff - 1] -= interp_quad_2(interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[n3_coeff - 2], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 2], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 2]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[n3_coeff - 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff - 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff - 1]),
																		 interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]));	
																		 // 486 = q2(483, 485, 487) = q2(c(183, 383, 583, 783), c(185, 385, 585, 785), c(187, 387, 587, 787))
					coeff_nodal_nodal_pos[n3_coeff] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff], nodal_nodal_nodal_pos[n3_coeff], nodal_nodal_nodal_pos[dim0_stride + n3_coeff], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff]);	// 487 = c(187, 387, 587, 787)
					if(even_n3){
						coeff_nodal_nodal_pos[n3_coeff + 1] -= interp_cubic(nodal_nodal_nodal_pos[-dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[n3_coeff + 1], nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1], nodal_nodal_nodal_pos[2*dim0_stride + n3_coeff + 1]);	// 488 = c(188, 388, 588, 788)
					}
				}
			}
			nodal_pos += dim0_stride;
            coeff_pos += dim0_stride;
		}
	}
	// decompse n1 x n2 x n3 data into coarse level (n1/2 x n2/2 x n3/2)
	void decompose_level_3D(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t dim0_stride, size_t dim1_stride){
		data_reorder_3D(data_pos, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
		compute_interpolant_difference_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n3_nodal = (n3 >> 1) + 1;
        compute_correction_3D(data_pos, data_buffer, load_v_buffer, n1, n2, n3, n1_nodal, h, dim0_stride, dim1_stride, default_batch_size);
        T * nodal_pos = data_pos;
        const T * correction_pos = data_buffer;
        for(int i=0; i<n1_nodal; i++){
            apply_correction_batched(nodal_pos, correction_pos, n2_nodal, dim1_stride, n3_nodal, true);
            nodal_pos += dim0_stride;
            correction_pos += n2_nodal * n3_nodal;
        }
	}
    void decompose_level_3D_with_hierarchical_basis(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t dim0_stride, size_t dim1_stride){
        data_reorder_3D(data_pos, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
        compute_interpolant_difference_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
    }
	void decompose_level_3D_cubic(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t dim0_stride, size_t dim1_stride){
        data_reorder_3D(data_pos, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
        compute_interpolant_difference_3D_cubic(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
    }
};

}

#endif