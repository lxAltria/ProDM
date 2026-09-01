#ifndef _MGARD_MISC_HPP
#define _MGARD_MISC_HPP

#include <string>

namespace MGARD{

template <class T>
void restriction_vertical_in_batch(T * data, size_t n1_nodal, size_t n1_coeff, size_t stride, int batchsize){
    T * nodal_pos = data;
    T * coeff_pos = data + n1_nodal * stride;
    // first nodal value
    for(int j=0; j<batchsize; j++){
        nodal_pos[j] += 0.5 * coeff_pos[j];
    }
    nodal_pos += stride;
    for(int i=1; i<n1_coeff; i++){
        for(int j=0; j<batchsize; j++){
            nodal_pos[j] += 0.5 * coeff_pos[j] + 0.5 * coeff_pos[j + stride];
        }
        nodal_pos += stride, coeff_pos += stride;
    }
    // last nodal value
    for(int j=0; j<batchsize; j++){
        nodal_pos[j] += 0.5 * coeff_pos[j];
    }
    // if n is even, the last nodal value can be skipped
}

template <class T>
void restriction_vertical(T * data, size_t n1, size_t n2, size_t stride, int batchsize=32){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    int num_batches = (n2_nodal - 1) / batchsize;
    T * data_pos = data;
    for(int i=0; i<num_batches; i++){
        restriction_vertical_in_batch(data_pos, n1_nodal, n1_coeff, stride, batchsize);
        data_pos += batchsize;
    }
    if(n2_nodal - batchsize * num_batches > 0){
        batchsize = n2_nodal - batchsize * num_batches;
        restriction_vertical_in_batch(data_pos, n1_nodal, n1_coeff, stride, batchsize);
    }
}

template <class T>
void restriction_1D(T * data, size_t n1){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    // backward pass
    T * nodal_pos = data;
    const T * coeff_pos = data + n1_nodal;
    nodal_pos[0] += 0.5 * coeff_pos[0];
    nodal_pos ++;
    for(int i=1; i<n1_coeff; i++){
        nodal_pos[0] += 0.5 * coeff_pos[0] + 0.5 * coeff_pos[1];
        nodal_pos ++, coeff_pos ++;
    }
    nodal_pos[0] += 0.5 * coeff_pos[0];
    // if n is even, the last nodal value can be skipped
}

template <class T>
void restriction_2D(T * data, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    // compute horizontal correction
    T * nodal_pos = data;
    for(int i=0; i<n1_nodal; i++){
        restriction_1D(nodal_pos, n2);
        nodal_pos += stride;
    }
    // compute vertical restriction
    restriction_vertical(data, n1, n2, stride);
}

template <class T>
void restriction_3D(T * data, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n2_nodal = (n2 >> 1) + 1;
    T * data_pos = data;
    // perform 2D restriction
    for(int i=0; i<n1; i++){
        restriction_2D(data_pos, n2, n3, dim1_stride);
        data_pos += dim0_stride;
    }        
    // compute vertical restriction
    data_pos = data;
    for(int i=0; i<n2_nodal; i++){
        restriction_vertical(data_pos, n1, n3, dim0_stride);
        data_pos += dim1_stride;
    }
}

template <class T>
void data_copy_3D(T * dst, const T * src, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    T * x_dst_pos = dst;
    const T * x_src_pos = src;
    for(int i=0; i<n1; i++){
        T * y_dst_pos = x_dst_pos;
        const T * y_src_pos = x_src_pos;
        for(int j=0; j<n2; j++){
            memcpy(y_dst_pos, y_src_pos, n3 * sizeof(T));
            y_dst_pos += dim1_stride;
            y_src_pos += dim1_stride;
        }
        x_dst_pos += dim0_stride;
        x_src_pos += dim1_stride;
    }
}

template <class T>
inline double dot_product(const T * A, const T * B, size_t n){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum += A[i] * B[i];
    }
    return sum;
}

template <class T>
double dot_product_3D(const T * A, const T * B, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    double sum = 0;
    const T * x_A_pos = A;
    const T * x_B_pos = B;
    for(int i=0; i<n1; i++){
        const T * y_A_pos = x_A_pos;
        const T * y_B_pos = x_B_pos;
        for(int j=0; j<n2; j++){
            // TODO: call CBLAS
            sum += dot_product(y_A_pos, y_B_pos, n3);
            y_A_pos += dim1_stride;
            y_B_pos += dim1_stride;
        }
        x_A_pos += dim0_stride;
        x_B_pos += dim1_stride;
    }
    return sum;
}

}
#endif