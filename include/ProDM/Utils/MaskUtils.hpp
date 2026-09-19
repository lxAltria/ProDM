#ifndef PRODM_UTILS_MASKUTILS_HPP
#define PRODM_UTILS_MASKUTILS_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "ProDM/Utils/RefactorUtils.hpp"
#include "ProDM/Compressor/ZSTD.hpp"

#include "ProDM/Namespace.hpp"

namespace ProDM {

template <class T>
void writemask(const char *filepath, T *data, size_t num_elements) {
    unsigned int bit_count = 1;
    unsigned int byte_count = bit_count / 8;
    unsigned int remainder_bit = bit_count % 8;
    size_t byteLength = 0;
    if (remainder_bit == 0) {
        byteLength = byte_count * num_elements + 1;
    }
    else {
        size_t tmp = remainder_bit * num_elements;
        byteLength = byte_count * num_elements + (tmp - 1) / 8 + 1;
    }
    std::vector<unsigned int> int_mask(num_elements, 0);
    for(int i=0; i<num_elements; i++){
        int_mask[i] = data[i];
    }
    std::vector<unsigned char> compressed_mask(byteLength, 0);
    if (byteLength != ProDM::save_fixed_length_bits(int_mask.data(), num_elements, compressed_mask.data(), bit_count)){}
    uint8_t * ZSTD_mask = nullptr;
    uint32_t ZSTD_mask_size = ZSTD::compress(compressed_mask.data(), byteLength, &ZSTD_mask);
    uint32_t mask_size = sizeof(size_t) + sizeof(uint32_t) + ZSTD_mask_size;
    uint8_t * mask_data = (uint8_t *) malloc(mask_size);
    uint8_t * mask_data_pos = mask_data;
    memcpy(mask_data_pos, &num_elements, sizeof(size_t));
    mask_data_pos += sizeof(size_t);
    memcpy(mask_data_pos, &ZSTD_mask_size, sizeof(uint32_t));
    mask_data_pos += sizeof(uint32_t);
    memcpy(mask_data_pos, ZSTD_mask, ZSTD_mask_size);
    FILE * file = fopen(filepath, "wb");
    if (file == nullptr) {
        perror("Error opening file");
        return;
    }
    fwrite(mask_data, 1, mask_size, file);
    fclose(file);
    free(mask_data);
    free(ZSTD_mask);
}

inline std::vector<unsigned char> readmask(const char *filepath, uint32_t & mask_file_size){
    FILE * file = fopen(filepath, "rb");
    if (file == nullptr){
        perror("Error opening file\n");
        return {};
    }
    fseek(file, 0, SEEK_END);
    uint32_t num_bytes = ftell(file);
    mask_file_size = num_bytes;
    rewind(file);
    uint8_t * mask_data = (uint8_t *)malloc(num_bytes);
    fread(mask_data, 1, num_bytes, file);
    fclose(file);
    uint8_t * mask_data_pos = mask_data;
    size_t num_elements;
    uint32_t ZSTD_mask_size;
    memcpy(&num_elements, mask_data_pos, sizeof(size_t));
    mask_data_pos += sizeof(size_t);
    memcpy(&ZSTD_mask_size, mask_data_pos, sizeof(uint32_t));
    mask_data_pos += sizeof(uint32_t);
    uint8_t * ZSTD_mask = (uint8_t *)malloc(ZSTD_mask_size);
    memcpy(ZSTD_mask, mask_data_pos, ZSTD_mask_size);
    mask_data_pos += ZSTD_mask_size;
    free(mask_data);
    size_t byteLength = 0;
    unsigned int bit_count = 1;
    unsigned int byte_count = bit_count / 8;
    unsigned remainder_bit = bit_count % 8;
    if (remainder_bit == 0){
        byteLength = byte_count * num_elements + 1;
    }
    else{
        size_t tmp = remainder_bit * num_elements;
        byteLength = byte_count * num_elements + (tmp - 1) / 8 + 1;
    }
    std::vector<unsigned char> compressed_mask(byteLength, 0);
    uint8_t * tmp_data = nullptr;
    uint32_t tmp_size = ZSTD::decompress(ZSTD_mask, ZSTD_mask_size, &tmp_data);

    if (tmp_data != nullptr){
        compressed_mask.assign(tmp_data, tmp_data + tmp_size);
        free(tmp_data);
    }
    std::vector<unsigned int> int_mask(num_elements, 0);
    std::vector<unsigned char> mask(num_elements, 0);
    if (compressed_mask.size() != ProDM::extract_fixed_length_bits(compressed_mask.data(), num_elements, int_mask.data(), bit_count)){}
    for(int i=0; i<num_elements; i++){
        mask[i] = int_mask[i];
    }
    free(ZSTD_mask);
    return mask;
}

}
#endif
