#ifndef _MDR_MGARD_DECOMPOSER_HPP
#define _MDR_MGARD_DECOMPOSER_HPP

#include "DecomposerInterface.hpp"
#include "decompose.hpp"
#include "recompose.hpp"
#include "decompose_new.hpp"
#include "recompose_new.hpp"
#include "decompose_interleave.hpp"
#include "reposition_recompose.hpp"
#include "decompose_interleave_new.hpp"
#include "reposition_recompose_new.hpp"
#include "coeff_decompose_interleave.hpp"
#include "coeff_reposition_recompose.hpp"

namespace MDR {
    // MGARD decomposer with orthogonal basis
    template<class T>
    class MGARDOrthoganalDecomposer : public concepts::DecomposerInterface<T> {
    public:
        MGARDOrthoganalDecomposer(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Decomposer<T> decomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                decomposer.decompose(data, dims, target_level, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                decomposer.decompose(data, dims, target_level, false, strs);
            }
        }
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Recomposer<T> recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                recomposer.recompose(data, dims, target_level, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                recomposer.recompose(data, dims, target_level, false, strs);
            }
        }
        void print() const {
            std::cout << "MGARD orthogonal decomposer" << std::endl;
        }
    };
    // MGARD decomposer with hierarchical basis
    template<class T>
    class MGARDHierarchicalDecomposer : public concepts::DecomposerInterface<T> {
    public:
        MGARDHierarchicalDecomposer(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Decomposer<T> decomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                decomposer.decompose(data, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                decomposer.decompose(data, dims, target_level, true, false, strs);
            }
        }
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Recomposer<T> recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                recomposer.recompose(data, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                recomposer.recompose(data, dims, target_level, true, false, strs);
            }
        }
        void print() const {
            std::cout << "MGARD hierarchical decomposer" << std::endl;
        }
    };
    // MGARD decomposer with cubic inperpolation
    template<class T>
    class MGARDCubicDecomposer : public concepts::DecomposerInterface<T>{
    public:
        MGARDCubicDecomposer(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Decomposer<T> decomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                decomposer.decompose(data, dims, target_level, false, true);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                decomposer.decompose(data, dims, target_level, false, true, strs);
            }
        }
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Recomposer<T> recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                recomposer.recompose(data, dims, target_level, false, true);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                recomposer.recompose(data, dims, target_level, false, true, strs);
            }
        }
        void print() const {
            std::cout << "MGARD cubic decomposer" << std::endl;
        }
    };
    // MGARD new decomposer with hierarchical basis
    template<class T>
    class MGARDHierarchicalDecomposer_new : public concepts::DecomposerInterface<T> {
    public:
        MGARDHierarchicalDecomposer_new(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Decomposer_new<T> decomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                decomposer.decompose(data, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                decomposer.decompose(data, dims, target_level, true, false, strs);
            }
        }
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {
            MGARD::Recomposer_new<T> recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                recomposer.recompose(data, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                recomposer.recompose(data, dims, target_level, true, false, strs);
            }
        }
        void print() const {
            std::cout << "MGARD hierarchical new decomposer" << std::endl;
        }
    };
    template<class T>
    class MGARDHierarchicalDecomposer_Interleaver : public concepts::DecomposerInterface<T> {
    public:
        MGARDHierarchicalDecomposer_Interleaver(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        // PSZ levels
        std::vector<std::vector<T>> decompose_interleave(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<std::vector<T>> level_buffers;
            MGARD::Decomposer_Interleaver<T> decomposer_interleaver;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                level_buffers = decomposer_interleaver.decompose(data, dims, target_level, true, false);
                level_buffer_dims = decomposer_interleaver.get_level_buffer_dims();
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                level_buffers = decomposer_interleaver.decompose(data, dims, target_level, true, false, strs);
                level_buffer_dims = decomposer_interleaver.get_level_buffer_dims();
            }
            return level_buffers;
        }
        // PSZ levels
        std::vector<T> reposition_recompose(std::vector<std::vector<T>>& level_buffers, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<T> data;
            MGARD::Repositioner_Recomposer_new<T> repositioner_recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false, strs);
            }
            return data;
        }
        // MGARD levels
        std::vector<std::vector<T>> decompose_interleave_combine_levels(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<std::vector<T>> level_buffers(target_level + 1);
            std::vector<std::vector<T>> level_buffers_;
            MGARD::Decomposer_Interleaver<T> decomposer_interleaver;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                level_buffers_ = decomposer_interleaver.decompose(data, dims, target_level, true, false);
                level_buffer_dims = decomposer_interleaver.get_level_buffer_dims();
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                level_buffers_ = decomposer_interleaver.decompose(data, dims, target_level, true, false, strs);
                level_buffer_dims = decomposer_interleaver.get_level_buffer_dims();
            }

            size_t num_dims = dims.size();
            level_buffers[0] = level_buffers_[0];
            for(int i=1; i<=target_level; i++){
                size_t level_x_size = 0;
                for(int j=1; j<=num_dims; j++){
                    level_x_size += level_buffers_[(i - 1) * num_dims + j].size();
                }
                level_buffers[i].resize(level_x_size);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
                T * level_x_buffers_pos = level_buffers[i].data();
                for(int j=1; j<=num_dims; j++){
                    size_t level_buffers_index = (i - 1) * num_dims + j;
                    // std::cout << "level_buffers_[" << level_buffers_index << "].size() = " << level_buffers_[level_buffers_index].size() << std::endl;
                    memcpy(level_x_buffers_pos, level_buffers_[level_buffers_index].data(), level_buffers_[level_buffers_index].size() * sizeof(T));
                    level_x_buffers_pos += level_buffers_[level_buffers_index].size();
                }
            }

            return level_buffers;
        }
        // MGARD levels
        std::vector<T> reposition_recompose_split_levels(std::vector<std::vector<T>>& level_buffers_, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<T> data;
            MGARD::Repositioner_Recomposer<T> repositioner_recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());

            size_t num_dims = dims.size();
            std::vector<std::vector<T>> level_buffers(target_level * num_dims + 1);
            level_buffers[0] = level_buffers_[0];
            auto level_dims = MGARD::compute_level_dims_new(dimensions, target_level);
            auto level_buffer_sizes = MGARD::compute_level_buffers_size(level_dims, target_level, level_buffer_dims);
            for(int i=1; i<level_buffer_sizes.size(); i++){
                level_buffers[i].resize(level_buffer_sizes[i]);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
            }
            for(int i=1; i<=target_level; i++){
                T * level_x_buffers_pos = level_buffers_[i].data();
                for(int j=1; j<=num_dims; j++){
                    size_t level_buffers_index = (i - 1) * num_dims + j;
                    // std::cout << "level_buffers_[" << level_buffers_index << "].size() = " << level_buffers_[level_buffers_index].size() << std::endl;
                    // std::cout << "level_buffer_sizes[" << level_buffers_index << "] = " << level_buffer_sizes[level_buffers_index] << std::endl;
                    memcpy(level_buffers[level_buffers_index].data(), level_x_buffers_pos, level_buffer_sizes[level_buffers_index] * sizeof(T));
                    level_x_buffers_pos += level_buffer_sizes[level_buffers_index];
                }
            }
            if(strides.size() == 0){
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false, strs);
            }
            return data;
        }
        void print() const {
            std::cout << "MGARD hierarchical decomposer & interleaver" << std::endl;
        }
        std::vector<std::vector<uint32_t>> get_level_buffer_dims(){
            return level_buffer_dims;
        }
    private:
        std::vector<std::vector<uint32_t>> level_buffer_dims;
    };
    template<class T>
    class MGARDHierarchicalDecomposer_Interleaver_new : public concepts::DecomposerInterface<T> {
    public:
        MGARDHierarchicalDecomposer_Interleaver_new(){}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        std::vector<std::vector<T>> decompose_interleave(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<std::vector<T>> level_buffers;
            MGARD::Decomposer_Interleaver_new<T> decomposer_interleaver;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                level_buffers = decomposer_interleaver.decompose(data, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                level_buffers = decomposer_interleaver.decompose(data, dims, target_level, true, false, strs);
            }
            return level_buffers;
        }
        std::vector<T> reposition_recompose(std::vector<std::vector<T>>& level_buffers, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<T> data;
            MGARD::Repositioner_Recomposer_new<T> repositioner_recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                data = repositioner_recomposer.recompose(level_buffers, dims, target_level, true, false, strs);
            }
            return data;
        }
        void print() const {
            std::cout << "MGARD hierarchical new decomposer & interleaver" << std::endl;
        }
    };
    template<class T>
    class MGARDHierarchical_Coeff_Decomposer_Interleaver : public concepts::DecomposerInterface<T> {
    public:
        MGARDHierarchical_Coeff_Decomposer_Interleaver(size_t direction_=0){
            direction = direction_;
        }
        MGARDHierarchical_Coeff_Decomposer_Interleaver() : MGARDHierarchical_Coeff_Decomposer_Interleaver(1) {}
        void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) const {}
        // PSZ levels
        std::vector<std::vector<T>> decompose_interleave(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<std::vector<T>> level_buffers;
            MGARD::Coeff_Decomposer_Interleaver<T> coeff_decomposer_interleaver;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                level_buffers = coeff_decomposer_interleaver.decompose(data, dims, direction, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                level_buffers = coeff_decomposer_interleaver.decompose(data, dims, direction, target_level, true, false, strs);
            }
            return level_buffers;
        }
        std::vector<T> reposition_recompose(std::vector<std::vector<T>>& level_buffers, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<T> data;
            MGARD::Coeff_Repositioner_Recomposer<T> coeff_repositioner_recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                data = coeff_repositioner_recomposer.recompose(level_buffers, dims, direction, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                data = coeff_repositioner_recomposer.recompose(level_buffers, dims, direction, target_level, true, false, strs);
            }
            return data;
        }
        // MGARD levels
        std::vector<std::vector<T>> decompose_interleave_combine_levels(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<std::vector<T>> level_buffers(target_level + 1);
            std::vector<std::vector<T>> level_buffers_;
            MGARD::Coeff_Decomposer_Interleaver<T> coeff_decomposer_interleaver;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());
            if(strides.size() == 0){
                level_buffers_ = coeff_decomposer_interleaver.decompose(data, dims, direction, target_level, true, false);
                level_buffer_dims = coeff_decomposer_interleaver.get_level_buffer_dims();
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                level_buffers_ = coeff_decomposer_interleaver.decompose(data, dims, direction, target_level, true, false, strs);
                level_buffer_dims = coeff_decomposer_interleaver.get_level_buffer_dims();
            }

            size_t num_dims = dims.size() - 1;
            level_buffers[0] = level_buffers_[0];
            for(int i=1; i<=target_level; i++){
                size_t level_x_size = 0;
                for(int j=1; j<=num_dims; j++){
                    level_x_size += level_buffers_[(i - 1) * num_dims + j].size();
                }
                level_buffers[i].resize(level_x_size);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
                T * level_x_buffers_pos = level_buffers[i].data();
                for(int j=1; j<=num_dims; j++){
                    size_t level_buffers_index = (i - 1) * num_dims + j;
                    // std::cout << "level_buffers_[" << level_buffers_index << "].size() = " << level_buffers_[level_buffers_index].size() << std::endl;
                    memcpy(level_x_buffers_pos, level_buffers_[level_buffers_index].data(), level_buffers_[level_buffers_index].size() * sizeof(T));
                    level_x_buffers_pos += level_buffers_[level_buffers_index].size();
                }
            }

            return level_buffers;
        }
        // MGARD levels
        std::vector<T> reposition_recompose_split_levels(std::vector<std::vector<T>>& level_buffers_, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides=std::vector<uint32_t>()) {
            std::vector<T> data;
            MGARD::Coeff_Repositioner_Recomposer<T> coeff_repositioner_recomposer;
            std::vector<size_t> dims(dimensions.begin(), dimensions.end());

            size_t num_dims = dims.size() - 1;
            std::vector<std::vector<T>> level_buffers(target_level * num_dims + 1);
            level_buffers[0] = level_buffers_[0];
            std::vector<uint32_t> dims_uint32;
            for(int i=0; i<dimensions.size(); i++){
                if(i != direction) dims_uint32.push_back(dimensions[i]);
            }
            auto level_dims = MGARD::compute_level_dims_new(dims_uint32, target_level);
            auto level_buffer_sizes = MGARD::compute_level_buffers_size(level_dims, target_level, level_buffer_dims);
            for(int i=0; i<level_buffer_sizes.size(); i++){
                level_buffer_sizes[i] *= dimensions[direction];
            }
            for(int i=1; i<level_buffer_sizes.size(); i++){
                level_buffers[i].resize(level_buffer_sizes[i]);
                // std::cout << "level_buffers[" << i << "].size() = " << level_buffers[i].size() << std::endl;
            }
            for(int i=1; i<=target_level; i++){
                T * level_x_buffers_pos = level_buffers_[i].data();
                for(int j=1; j<=num_dims; j++){
                    size_t level_buffers_index = (i - 1) * num_dims + j;
                    // std::cout << "level_buffers_[" << level_buffers_index << "].size() = " << level_buffers_[level_buffers_index].size() << std::endl;
                    // std::cout << "level_buffer_sizes[" << level_buffers_index << "] = " << level_buffer_sizes[level_buffers_index] << std::endl;
                    memcpy(level_buffers[level_buffers_index].data(), level_x_buffers_pos, level_buffer_sizes[level_buffers_index] * sizeof(T));
                    level_x_buffers_pos += level_buffer_sizes[level_buffers_index];
                }
            }
            if(strides.size() == 0){
                data = coeff_repositioner_recomposer.recompose(level_buffers, dims, direction, target_level, true, false);
            }
            else{
                std::vector<size_t> strs(strides.begin(), strides.end());
                data = coeff_repositioner_recomposer.recompose(level_buffers, dims, direction, target_level, true, false, strs);
            }
            return data;
        }
        void print() const {
            std::cout << "MGARD hierarchical new decomposer & interleaver" << std::endl;
        }
    private:
        size_t direction;
        std::vector<std::vector<uint32_t>> level_buffer_dims;
    };
}
#endif
