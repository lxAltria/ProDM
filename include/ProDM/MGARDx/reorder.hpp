#ifndef _MGARD_REORDER_HPP
#define _MGARD_REORDER_HPP

#include <vector>

namespace MGARD{

using namespace std;

// switch the rows in the data for coherent memory access
/*
@params data_pos: starting position of data
@params data_buffer: buffer to store intermediate data
@params n1, n2: dimensions
@params stride: stride for the non-continguous dimension
Illustration:
o: nodal data, x: coefficient data
    oooxx       oooxx
    xxxxx       oooxx
    oooxx   =>  oooxx
    xxxxx       xxxxx
    oooxx       xxxxx
*/
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

// inverse operation for switch_rows_2D_by_buffer
/*
@params data_pos: starting position of data
@params data_buffer: buffer to store intermediate data
@params n1, n2: dimensions
@params stride: stride for the non-continguous dimension
Illustration:
o: nodal data, x: coefficient data
    oooxx       oooxx
    oooxx       xxxxx
    oooxx   =>  oooxx
    xxxxx       xxxxx
    xxxxx       oooxx
*/
template <class T>
void switch_rows_2D_by_buffer_reverse(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    T * nodal_data_buffer = data_pos + stride; // skip the first nodal row
    T * coeff_data_buffer = data_pos + n1_nodal * stride;
    T * cur_data_pos = data_buffer + n2;
    for(int i=0; i<n1_coeff; i++){
        // copy coefficient rows
        memcpy(cur_data_pos, coeff_data_buffer + i * stride, n2 * sizeof(T));
        cur_data_pos += n2;
        // copy nodal rows
        memcpy(cur_data_pos, nodal_data_buffer + i * stride, n2 * sizeof(T));
        cur_data_pos += n2;
    }
    if(!(n1&1)){
        // n1 is even, move the last nodal row
        memcpy(cur_data_pos, coeff_data_buffer - stride, n2 * sizeof(T));
    }
    // copy data back
    cur_data_pos = data_pos + stride;
    for(int i=1; i<n1; i++){
        memcpy(cur_data_pos, data_buffer + i * n2, n2 * sizeof(T));
        cur_data_pos += stride;
    }
}

// reorder the data to put all the coefficient to the back
template <class T>
void data_reorder_1D(const T * data_pos, size_t n_nodal, size_t n_coeff, T * nodal_buffer, T * coeff_buffer){
    T * nodal_pos = nodal_buffer;
    T * coeff_pos = coeff_buffer;
    T const * cur_data_pos = data_pos;
    for(int i=0; i<n_coeff; i++){
        *(nodal_pos++) = *(cur_data_pos++);
        *(coeff_pos++) = *(cur_data_pos++);
    }
    *(nodal_pos++) = *(cur_data_pos++);
    if(n_nodal == n_coeff + 2){
        // if even, add a nodal value such that the interpolant
        // of the last two nodal values equal to the last coefficient
        *nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
        // *nodal_pos -= nodal_pos[-1];
    }
}

// reorder the data to put all the coefficient to the back
/*
    oxoxo       oooxx       oooxx
    xxxxx   (1) xxxxx   (2) oooxx
    oxoxo   =>  oooxx   =>  oooxx
    xxxxx       xxxxx       xxxxx
    oxoxo       oooxx       xxxxx
*/
template <class T>
void data_reorder_2D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    for(int i=0; i<n1; i++){
        data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff row into nodal row
        cur_data_pos -= stride;
        for(int j=0; j<n2; j++){
            cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[-stride + j];
            // cur_data_pos[j] -= cur_data_pos[-stride + j];
        }
    }
    // do reorder (2)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
}

template <class T>
void data_reorder_2D_new(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    for(int i=0; i<n1; i+=2){
        data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        cur_data_pos += 2*stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff row into nodal row
        cur_data_pos -= stride;
        data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        for(int j=0; j<n2; j++){
            cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[-stride + j];
        }
    }
    // do reorder (2)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
}

/*
    2D reorder + vertical reorder
*/
template <class T>
void data_reorder_3D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // do 2D reorder
    for(int i=0; i<n1; i++){
        data_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff plane into nodal plane
        cur_data_pos -= dim0_stride;
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                cur_data_pos[k] = 2 * cur_data_pos[k] - cur_data_pos[- dim0_stride + k];
                // cur_data_pos[k] -= cur_data_pos[- dim0_stride + k];
            }
            cur_data_pos += dim1_stride;
        }
    }
    cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
}

template <class T>
void data_reorder_3D_new(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // do 2D reorder
    for(int i=0; i<n1; i+=2){
        data_reorder_2D_new(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += 2*dim0_stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff plane into nodal plane
        cur_data_pos -= dim0_stride;
        data_reorder_2D_new(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                cur_data_pos[k] = 2 * cur_data_pos[k] - cur_data_pos[- dim0_stride + k];
            }
            cur_data_pos += dim1_stride;
        }
    }
    cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
}

// reorder the data to original order (insert coeffcients between nodal values)
template <class T>
void data_reverse_reorder_1D(T * data_pos, int n_nodal, int n_coeff, const T * nodal_buffer, const T * coeff_buffer){
    const T * nodal_pos = nodal_buffer;
    const T * coeff_pos = coeff_buffer;
    T * cur_data_pos = data_pos;
    for(int i=0; i<n_coeff; i++){
        *(cur_data_pos++) = *(nodal_pos++);
        *(cur_data_pos++) = *(coeff_pos++);
    }
    *(cur_data_pos++) = *(nodal_pos++);
    if(n_nodal == n_coeff + 2){
        // if even, the last coefficient equals to the interpolant
        // of the last two nodal values
        *cur_data_pos = (nodal_pos[-1] + nodal_pos[0]) / 2;
    }
}

template <class T>
void data_reverse_reorder_1D_for_generating_purpose(T * data_pos, int n_nodal, int n_coeff, const T * nodal_buffer, const T * coeff_buffer){
    const T * nodal_pos = nodal_buffer;
    const T * coeff_pos = coeff_buffer;
    T * cur_data_pos = data_pos;
    for(int i=0; i<n_coeff; i++){
        *(cur_data_pos++) = *(nodal_pos++);
        *(cur_data_pos++) = *(coeff_pos++);
    }
    *(cur_data_pos++) = *(nodal_pos++);
}

/*
    oooxx       oooxx       oxoxo
    oooxx   (1) xxxxx   (2) xxxxo
    oooxx   =>  oooxx   =>  oxoxo
    xxxxx       xxxxx       xxxxx
    xxxxx       oooxx       oxoxo
*/
template <class T>
void data_reverse_reorder_2D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer_reverse(data_pos, data_buffer, n1, n2, stride);
    // do reorder (2)
    for(int i=0; i<n1; i++){
        memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
        data_reverse_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        cur_data_pos += stride;
    }
    if(!(n1 & 1)){
        // n1 is even, recover the coefficients
        cur_data_pos -= stride;
        for(int j=0; j<n2; j++){
            cur_data_pos[j] = (cur_data_pos[j] + cur_data_pos[-stride + j]) / 2;
        }
    }
}

template <class T>
void data_reverse_reorder_2D_for_generating_purpose(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer_reverse(data_pos, data_buffer, n1, n2, stride);
    // do reorder (2)
    for(int i=0; i<n1; i++){
        memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
        data_reverse_reorder_1D_for_generating_purpose(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        cur_data_pos += stride;
    }
}

template <class T>
void data_reverse_reorder_2D_new(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    // TODO: change to online processing for memory saving
    switch_rows_2D_by_buffer_reverse(data_pos, data_buffer, n1, n2, stride);
    // do reorder (2)
    for(int i=0; i<n1; i+=2){
        memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
        data_reverse_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        cur_data_pos += 2*stride;
    }
    if(!(n1 & 1)){
        // n1 is even, recover the coefficients
        cur_data_pos -= stride;
        memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
        data_reverse_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        for(int j=0; j<n2; j++){
            cur_data_pos[j] = (cur_data_pos[j] + cur_data_pos[-stride + j]) / 2;
        }
    }
}

/*
    vertical reorder + 2D reorder
*/
template <class T>
void data_reverse_reorder_3D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer_reverse(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
    // do 2D reorder
    cur_data_pos = data_pos;
    for(int i=0; i<n1; i++){
        data_reverse_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff plane into nodal plane
        cur_data_pos -= dim0_stride;
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                cur_data_pos[k] = (cur_data_pos[k] + cur_data_pos[- dim0_stride + k]) / 2;
            }
            cur_data_pos += dim1_stride;
        }
    }
}

template <class T>
void data_reverse_reorder_3D_for_generating_purpose(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer_reverse(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
    // do 2D reorder
    cur_data_pos = data_pos;
    for(int i=0; i<n1; i++){
        data_reverse_reorder_2D_for_generating_purpose(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
}

template <class T>
void data_reverse_reorder_3D_new(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        switch_rows_2D_by_buffer_reverse(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
    // do 2D reorder
    cur_data_pos = data_pos;
    for(int i=0; i<n1; i+=2){
        data_reverse_reorder_2D_new(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += 2*dim0_stride;
    }
    if(!(n1 & 1)){
        // n1 is even, change the last coeff plane into nodal plane
        cur_data_pos -= dim0_stride;
        data_reverse_reorder_2D_new(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                cur_data_pos[k] = (cur_data_pos[k] + cur_data_pos[- dim0_stride + k]) / 2;
            }
            cur_data_pos += dim1_stride;
        }
    }
}

}
#endif