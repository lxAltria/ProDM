#ifndef _MDR_REFACTOR_UTILS_HPP
#define _MDR_REFACTOR_UTILS_HPP

#include <cassert>
#include <vector>
#include <cmath>
#include <ctime>

namespace MDR {

    // MDR utility functions

    // MGARD related
    // TODO: put API in MGARD

    // compute level dimensions
    /*
        @params dims: input dimensions
        @params target_level: the target decomposition level
    */
    std::vector<std::vector<uint32_t>> compute_level_dims(const std::vector<uint32_t>& dims, uint32_t target_level){
        std::vector<std::vector<uint32_t>> level_dims;
        for(int i=0; i<=target_level; i++){
            level_dims.push_back(std::vector<uint32_t>(dims.size()));
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

    // compute level elements
    /*
        @params level_dims: dimensions for all levels
        @params target_level: the target decomposition level
    */
    std::vector<uint32_t> compute_level_elements(const std::vector<std::vector<uint32_t>>& level_dims, int target_level){
        assert(level_dims.size());
        uint8_t num_dims = level_dims[0].size();
        std::vector<uint32_t> level_elements(level_dims.size());
        level_elements[0] = 1;
        for(int j=0; j<num_dims; j++){
            level_elements[0] *= level_dims[0][j];
        }
        uint32_t pre_num_elements = level_elements[0];
        for(int i=1; i<=target_level; i++){
            uint32_t num_elements = 1;
            for(int j=0; j<num_dims; j++){
                num_elements *= level_dims[i][j];
            }
            level_elements[i] = num_elements - pre_num_elements;
            pre_num_elements = num_elements;
        }
        return level_elements;
    }

    std::vector<uint32_t> compute_level_buffers_size(const std::vector<std::vector<uint32_t>>& level_dims, int target_level){
        assert(level_dims.size());
        size_t num_dims = level_dims[0].size();
        // size_t count;
        size_t size;
        std::vector<uint32_t> level_sizes(1 + num_dims * target_level);
        level_sizes[0] = 1;
        // std::cout << "level 0, dims = ";
        for(int i=0; i<num_dims; i++){
            level_sizes[0] *= level_dims[0][i];
            // std::cout << level_dims[0][i] << " ";
        }
        // std::cout << std::endl;
        for(size_t l=1; l<=target_level; l++){
            // count = 0;
            switch (num_dims){
                case 1:
                {
                    size = level_dims[l][0];
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    // assert(count == (level_dims[l][0] - level_dims[l-1][0]));
                    break;
                }
                case 2:
                {
                    size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1];
                    // std::cout << "level " << (l-1) * num_dims + 1 << ", dims = " << level_dims[l][0] - level_dims[l-1][0] << " " << level_dims[l-1][1] << std::endl;
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]);
                    // std::cout << "level " << (l-1) * num_dims + 2 << ", dims = " << level_dims[l][0] << " " << level_dims[l][1] - level_dims[l-1][1] << std::endl;
                    level_sizes[(l-1) * num_dims + 2] = size;
                    // count += size;
                    // assert(count == (level_dims[l][0]*level_dims[l][1] - level_dims[l-1][0]*level_dims[l-1][1]));
                    break;
                }
                case 3:
                {
                    // n1 (cur_n1 - pre_n1) * pre_n2 * pre_n3
                    size = (level_dims[l][0] - level_dims[l-1][0]) * level_dims[l-1][1] * level_dims[l-1][2];
                    // std::cout << "level " << (l-1) * num_dims + 1 << ", dims = " << (level_dims[l][0] - level_dims[l-1][0]) << " " << level_dims[l-1][1] << " " << level_dims[l-1][2] << std::endl;
                    level_sizes[(l-1) * num_dims + 1] = size;
                    // count += size;
                    // n2 cur_n1 * (cur_n2 - pre_n2) * pre_n3
                    size = level_dims[l][0] * (level_dims[l][1] - level_dims[l-1][1]) * level_dims[l-1][2];
                    // std::cout << "level " << (l-1) * num_dims + 2 << ", dims = " << level_dims[l][0] << " " << (level_dims[l][1] - level_dims[l-1][1]) << " " << level_dims[l-1][2] << std::endl;
                    level_sizes[(l-1) * num_dims + 2] = size;
                    // count += size;
                    // n3 cur_n1 * cur_n2 * (cur_n3 - pre_n3)
                    size = level_dims[l][0] * level_dims[l][1] * (level_dims[l][2] - level_dims[l-1][2]);
                    // std::cout << "level " << (l-1) * num_dims + 3 << ", dims = " << level_dims[l][0] << " " << level_dims[l][1] << " " << (level_dims[l][2] - level_dims[l-1][2]) << std::endl;
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

    // Simple utility functions

    // compute maximum value in level
    /*
    @params data: level data
    @params n: number of level data points
    */
    template <class T>
    T compute_max_abs_value(const T * data, uint32_t n){
        T max_val = 0;
        for(int i=0; i<n; i++){
            T val = fabs(data[i]);
            if(val > max_val) max_val = val;
        }
        return max_val;
    }

    // Get size of vector
    template <class T>
    inline uint32_t get_size(const std::vector<T>& vec){
        return vec.size() * sizeof(T);
    }
    template <class T>
    uint32_t get_size(const std::vector<std::vector<T>>& vec){
        uint32_t size = 0;
        for(int i=0; i<vec.size(); i++){
            size += sizeof(uint32_t) + vec[i].size() * sizeof(T);
        }
        return size;
    }

    // Serialize/deserialize vectors
    // Auto-increment buffer position
    template <class T>
    inline void serialize(const std::vector<T>& vec, uint8_t *& buffer_pos){
        memcpy(buffer_pos, vec.data(), vec.size() * sizeof(T));
        buffer_pos += vec.size() * sizeof(T);
    }
    template <class T>
    void serialize(const std::vector<std::vector<T>>& vec, uint8_t *& buffer_pos){
        uint8_t const * const start = buffer_pos;
        for(int i=0; i<vec.size(); i++){
            *reinterpret_cast<uint32_t*>(buffer_pos) = vec[i].size();
            buffer_pos += sizeof(uint32_t);
            memcpy(buffer_pos, vec[i].data(), vec[i].size() * sizeof(T));
            buffer_pos += vec[i].size() * sizeof(T);
        }
    }
    template <class T>
    inline void deserialize(uint8_t const *& buffer_pos, uint32_t size, std::vector<T>& vec){
        vec.clear();
        vec = std::vector<T>(reinterpret_cast<const T*>(buffer_pos), reinterpret_cast<const T*>(buffer_pos) + size);
        buffer_pos += size * sizeof(T);
    }
    template <class T>
    void deserialize(uint8_t const *& buffer_pos, uint32_t num_levels, std::vector<std::vector<T>>& vec){
        vec.clear();
        for(int i=0; i<num_levels; i++){
            uint32_t num = *reinterpret_cast<const uint32_t*>(buffer_pos);
            buffer_pos += sizeof(uint32_t);
            std::vector<T> level_vec = std::vector<T>(reinterpret_cast<const T *>(buffer_pos), reinterpret_cast<const T *>(buffer_pos) + num);
            vec.push_back(level_vec);
            buffer_pos += num * sizeof(T);
        }
    }

    // print vector
    template <class T>
    void print_vec(const std::vector<T>& vec){
        for(int i=0; i<vec.size(); i++){
            std::cout << vec[i] << " ";
        }
        std::cout << std::endl;
    }
    // print nested vector
    template <class T>
    void print_vec(const std::string& name, const std::vector<std::vector<T>>& vec){
        std::cout << name << std::endl;
        for(int i=0; i<vec.size(); i++){
            print_vec(vec[i]);
        }
        std::cout << std::endl;
    }

    class Timer{
    public:
        void start(){
            err = clock_gettime(CLOCK_REALTIME, &start_time);
        }
        void end(){
            err = clock_gettime(CLOCK_REALTIME, &end_time);
            total_time += (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/(double)1000000000;
        }
        double get(){
            double time = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/(double)1000000000;
            clear();
            return time;
        }
        void clear(){
            total_time = 0;
        }
        void print(std::string s){
            std::cout << s << " time: " << total_time << "s" << std::endl;
            clear();
        }
    private:
        int err = 0;
        double total_time = 0;
        struct timespec start_time, end_time;
    };

}
#endif
