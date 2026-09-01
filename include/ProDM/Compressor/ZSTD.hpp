#ifndef _MDR_ZSTD_HPP
#define _MDR_ZSTD_HPP

#include "zstd.h"
#include <cstdlib>
#include <iostream>

namespace MDR {
    namespace ZSTD{
        #define ZSTD_LEVEL 3 //default setting of level is 3
        // ZSTD lossless compressor
        uint32_t compress(const uint8_t* data, uint32_t dataLength, uint8_t** compressBytes) {
            size_t dstCapacity = ZSTD_compressBound(dataLength);
            *compressBytes = (uint8_t*)malloc(dstCapacity + sizeof(size_t));
            *reinterpret_cast<size_t*>(*compressBytes) = dataLength;
            size_t outSize = ZSTD_compress(*compressBytes + sizeof(size_t), dstCapacity, data, dataLength, ZSTD_LEVEL);
            if(ZSTD_isError(outSize)){
                std::cerr << "ZSTD compression error: " << ZSTD_getErrorName(outSize) << std::endl;
                exit(-1);
            }
            return outSize + sizeof(size_t);
        }
        uint32_t decompress(const uint8_t* compressBytes, uint32_t cmpSize, uint8_t** oriData) {
            size_t outSize = *reinterpret_cast<const size_t*>(compressBytes);
            *oriData = (uint8_t*)malloc(outSize);
            size_t decompressedSize = ZSTD_decompress(*oriData, outSize, compressBytes + sizeof(size_t), cmpSize - sizeof(size_t));
            if(ZSTD_isError(decompressedSize)){
                std::cerr << "ZSTD decompression error: " << ZSTD_getErrorName(decompressedSize) << std::endl;
                exit(-1);
            }
            return outSize;
        }
    }
}
#endif
