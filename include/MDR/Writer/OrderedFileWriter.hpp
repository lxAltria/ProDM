#ifndef _MDR_SERIAL_FILE_WRITER_HPP
#define _MDR_SERIAL_FILE_WRITER_HPP

#include "WriterInterface.hpp"
#include <cstdio>

namespace MDR {
    // A writer that writes the serialized components
    class OrderedFileWriter : public concepts::WriterInterface {
    public:
        OrderedFileWriter(const std::string& metadata_file, const std::string& data_file) : metadata_file(metadata_file), data_file(data_file) {}

        std::vector<uint32_t> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes) const {
            std::vector<uint32_t> temp;
            return temp;
        }

        uint32_t write_components(uint8_t const * data, uint32_t size) const {
            FILE * file = fopen(data_file.c_str(), "wb");
            fwrite(data, 1, size, file);
            fclose(file);
            return 0;
        }

        void write_metadata(uint8_t const * metadata, uint32_t size) const {
            FILE * file = fopen(metadata_file.c_str(), "wb");
            fwrite(metadata, 1, size, file);
            fclose(file);
        }

        ~OrderedFileWriter(){}

        void print() const {
            std::cout << "Ordered file writer." << std::endl;
        }
    private:
        std::string metadata_file;
        std::string data_file;
    };
}
#endif
