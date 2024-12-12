#ifndef _MDR_MGARD_DECOMPOSER_HPP
#define _MDR_MGARD_DECOMPOSER_HPP

#include "DecomposerInterface.hpp"
#include "decompose.hpp"
#include "recompose.hpp"

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
}
#endif
