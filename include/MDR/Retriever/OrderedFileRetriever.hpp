#ifndef _MDR_ORDERED_FILE_RETRIEVER_HPP
#define _MDR_ORDERED_FILE_RETRIEVER_HPP

#include "RetrieverInterface.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>

namespace MDR {

    // Data retriever for files
    class OrderedFileRetriever : public concepts::RetrieverInterface {
    public:
        OrderedFileRetriever(const std::string& metadata_file,
                             const std::string& data_file)
            : metadata_file(metadata_file),
              data_file(data_file),
              offset(0),
              components(nullptr),
              total_retrieved_size(0) {}

        std::vector<std::vector<const uint8_t*>> retrieve_level_components(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<uint32_t>& retrieve_sizes, const std::vector<uint8_t>& prev_level_num_bitplanes, const std::vector<uint8_t>& level_num_bitplanes){
            std::vector<std::vector<const uint8_t*>> tmp;
            return tmp;
        }

        uint8_t* retrieve_components(uint32_t retrieve_size) {
            FILE* file = fopen(data_file.c_str(), "rb");
            if (!file) {
                std::cerr << "Failed to open data file: " << data_file << std::endl;
                return nullptr;
            }

            if (fseek(file, static_cast<long>(offset), SEEK_SET)) {
                std::cerr << "Errors in fseek while retrieving from file" << std::endl;
            }

            uint8_t* buffer = static_cast<uint8_t*>(std::malloc(retrieve_size));
            if (!buffer) {
                std::cerr << "Failed to allocate buffer in retrieve_components" << std::endl;
                fclose(file);
                return nullptr;
            }

            size_t read_bytes = std::fread(buffer, sizeof(uint8_t), retrieve_size, file);
            if (read_bytes != retrieve_size) {
                std::cerr << "Warning: requested " << retrieve_size
                          << " bytes but only read " << read_bytes << " bytes." << std::endl;
            }
            fclose(file);

            components = buffer;

            offset += static_cast<uint32_t>(read_bytes);
            total_retrieved_size += read_bytes;

            return buffer;
        }

        uint8_t* load_metadata() const {
            FILE* file = fopen(metadata_file.c_str(), "rb");
            if (!file) {
                std::cerr << "Failed to open metadata file: " << metadata_file << std::endl;
                return nullptr;
            }
            fseek(file, 0, SEEK_END);
            uint32_t num_bytes = static_cast<uint32_t>(ftell(file));
            rewind(file);

            uint8_t* metadata = static_cast<uint8_t*>(std::malloc(num_bytes));
            if (!metadata) {
                std::cerr << "Failed to allocate metadata buffer" << std::endl;
                fclose(file);
                return nullptr;
            }

            size_t read_bytes = std::fread(metadata, 1, num_bytes, file);
            if (read_bytes != num_bytes) {
                std::cerr << "Warning: requested " << num_bytes
                          << " bytes of metadata but only read " << read_bytes << " bytes."
                          << std::endl;
            }
            fclose(file);
            return metadata;
        }

        void release() {
            if (components != nullptr) {
                std::free(components);
                components = nullptr;
            }
        }

        size_t get_retrieved_size() {
            return total_retrieved_size;
        }

        uint32_t get_offset() {
            return offset;
        }

        ~OrderedFileRetriever() {
            release();
        }

        void print() const {
            std::cout << "Ordered file retriever." << std::endl;
        }

    private:
        std::string metadata_file;
        std::string data_file;
        uint32_t offset;
        uint8_t* components;
        size_t total_retrieved_size;
    };

} // namespace MDR

#endif
