#ifndef _MDR_NOMASK_GE_SYNTHESIZER_HPP
#define _MDR_NOMASK_GE_SYNTHESIZER_HPP

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <numeric>
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "MDR/Refactor/Refactor.hpp"
#include "SZ3/api/sz.hpp"
#include "PDR/Refactor/Refactor.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
// #include "PDR/Reconstructor/ApproximationBasedReconstructor.hpp"
// #include "PDR/Reconstructor/WeightedApproximationBasedReconstructor.hpp"
#include "ompSZp_typemanager.h"
#include <cstdint>

const std::vector<std::string> varlist = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
const int n_vars = 5;

namespace MDR {

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer> generateRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer, bool negabinary){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    return refactor;
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

template <class T, class Decomposer, class InterleaverT, class InterleaverInt, class Encoder, class Compressor, class ErrorCollector, class Writer>
MDR::QoIRefactor<T, Decomposer, InterleaverT, InterleaverInt, Encoder, Compressor, ErrorCollector, Writer> generateRefactor(Decomposer decomposer, InterleaverT interleaver, InterleaverInt weight_interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer, bool negabinary){
    auto refactor = MDR::QoIRefactor<T, Decomposer, InterleaverT, InterleaverInt, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, weight_interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    return refactor;
}

template <class T, class Approximator, class Encoder, class Compressor, class Writer>
PDR::WeightedApproximationBasedRefactor<T, Approximator, Encoder, Compressor, Writer> generateWBPRefactor(Approximator approximator, Encoder encoder, Compressor compressor, Writer writer, bool negabinary){
    auto refactor = PDR::WeightedApproximationBasedRefactor<T, Approximator, Encoder, Compressor, Writer>(approximator, encoder, compressor, writer);
    refactor.negabinary = negabinary;
    return refactor;
}

template <class T, class Approximator, class Encoder, class Compressor, class Writer>
PDR::ApproximationBasedRefactor<T, Approximator, Encoder, Compressor, Writer> generateBPRefactor(Approximator approximator, Encoder encoder, Compressor compressor, Writer writer, bool negabinary){
    auto refactor = PDR::ApproximationBasedRefactor<T, Approximator, Encoder, Compressor, Writer>(approximator, encoder, compressor, writer);
    refactor.negabinary = negabinary;
    return refactor;
}

template <class T, class Decomposer, class InterleaverT, class InterleaverInt, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
MDR::WeightReconstructor<T, Decomposer, InterleaverT, InterleaverInt, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateReconstructor(Decomposer decomposer, InterleaverT interleaver, InterleaverInt weight_interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::WeightReconstructor<T, Decomposer, InterleaverT, InterleaverInt, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, weight_interleaver, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

template <class T, class Approximator, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
PDR::WeightedApproximationBasedReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateWBPReconstructor(Approximator approximator, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = PDR::WeightedApproximationBasedReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(approximator, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

template <class T, class Approximator, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
PDR::ApproximationBasedReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateBPReconstructor(Approximator approximator, Encoder encoder, Compressor compressor,
 ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = PDR::ApproximationBasedReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(approximator, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

template <class T, class Approximator, class Encoder, class Compressor, class Writer>
PDR::GERefactor<T, Approximator, Encoder, Compressor, Writer> generateWBPRefactor_GE(Approximator approximator, Encoder encoder, Compressor compressor, Writer writer, bool negabinary){
    auto refactor = PDR::GERefactor<T, Approximator, Encoder, Compressor, Writer>(approximator, encoder, compressor, writer);
    refactor.negabinary = negabinary;
    return refactor;
}

template <class T, class Approximator, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
PDR::GEReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateWBPReconstructor_GE(Approximator approximator, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = PDR::GEReconstructor<T, Approximator, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(approximator, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

const int target_level = 8;
const int num_bitplanes = 60;

template <class T>
static T compute_vr(const std::vector<T>& vec){
    T min = vec[0];
    T max = vec[0];
    for(int i=0; i<vec.size(); i++){
        if(vec[i] < min) min = vec[i];
        if(vec[i] > max) max = vec[i];
    }
    return max - min;
}

template<class T>
char * SZ3_compress(size_t num_elements, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(num_elements);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
char * SZ3_compress_3D(size_t num_elements, uint32_t n1, uint32_t n2, uint32_t n3, T * data, double abs_eb, size_t& compressed_size){
    SZ3::Config conf(n1, n2, n3);
    conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs_eb;
    size_t cmpSize = 0;
    char *cmpData = SZ_compress<T>(conf, data, cmpSize);
    compressed_size = cmpSize;
    return cmpData;
}

template<class T>
void SZ3_decompress(char * cmpData, size_t compressed_size, T * dec_data){
    SZ3::Config conf1;
    SZ_decompress<T>(conf1, cmpData, compressed_size, dec_data);
}

inline int find_index(double target_rel_eb, double& rel_eb){
    int i = 0;
    while(target_rel_eb < rel_eb){
        i ++;
        rel_eb /= 10;
    }
    return i;
}

template<class Type>
void refactor_GE(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<Type>();
        auto interleaver = MDR::DirectInterleaver<Type>();
        auto encoder = MDR::PerBitBPEncoder<Type, uint32_t>();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<Type>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<Type>(decomposer, interleaver, encoder, compressor, collector, writer);
        if(i < 3){
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            std::cout << "index = " << index << std::endl;
            refactor.refactor(buffer.data(), dims_masked, target_level, num_bitplanes);            
        } 
        else{
            refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
        }
    }
}

template<class Type>
void refactor_GE_SZ3(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    std::vector<double> rel_ebs;
    const int num_snapshot = 18;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored/SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_valid_data, buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
            std::cout << "index = " << index << std::endl;
        } 
        else{
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored/SZ3_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_elements, vars_vec[i].data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                free(compressed_data);
            }
        }
    }
}

template<class Type>
void refactor_GE_SZ3_delta(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    // std::string mask_file = rdata_file_prefix + "mask.bin";
    // MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    std::vector<std::vector<Type>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<uint32_t> dims_masked;
    dims_masked.push_back(num_valid_data);
    std::vector<Type> buffer(num_valid_data);

    std::vector<double> value_range(n_vars);
    for(int i=0; i<n_vars; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    std::vector<double> rel_ebs;
    const int num_snapshot = 18;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_vars; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        if(i < 3){
            // use masked refactoring for vx vy vz
            int index = 0;
            for(int j=0; j<num_elements; j++){
                if(mask[j]){
                    buffer[index ++] = vars_vec[i][j];
                }
            }
            std::vector<Type> data_buffer(buffer);
            std::vector<Type> dec_data_buffer(buffer);
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_valid_data, data_buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
                for(int i=0; i<num_valid_data; i++){
                    data_buffer[i] = data_buffer[i] - dec_data_buffer[i];
                }
                free(compressed_data);
            }
            std::cout << "index = " << index << std::endl;
        } 
        else{
            std::vector<Type> data_buffer(vars_vec[i]);
            std::vector<Type> dec_data_buffer(vars_vec[i]);
            for(int j=0; j<num_snapshot; j++){
                std::string filename = rdir_prefix + "_refactored/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                size_t compressed_size = 0;
                auto compressed_data = SZ3_compress(num_elements, data_buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
                MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
                SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
                for(int i=0; i<num_elements; i++){
                    data_buffer[i] = data_buffer[i] - dec_data_buffer[i];
                }
                free(compressed_data);
            }
        }
    }
}

template<class Type>
void refactor_S3D(uint32_t n1, uint32_t n2, uint32_t n3, const std::string s3d_data_file_prefix, const std::string s3d_rdata_file_prefix){
    std::vector<std::string> species = {"H2", "O2", "H2O", "H", "O", "OH"};
    int n_species = species.size();

    std::vector<std::vector<Type>> vars_vec(n_species);

    size_t num_elements = 0;
    for(int i=0; i<n_species; i++){
        auto Xi = MGARD::readfile<Type>((s3d_data_file_prefix + species[i] + ".dat").c_str(), num_elements);
        vars_vec[i] = Xi;
    }

    int target_level = 4;
    std::vector<uint32_t> dims = {n1, n2, n3};

    for(int i=0; i<n_species; i++){
        std::string rdir_prefix = s3d_rdata_file_prefix + species[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<Type>();
        auto interleaver = MDR::DirectInterleaver<Type>();
        auto encoder = MDR::PerBitBPEncoder<Type, uint32_t>();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<Type>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<Type>(decomposer, interleaver, encoder, compressor, collector, writer);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);  
    }
}

template<class Type>
void refactor_S3D_SZ3(uint32_t n1, uint32_t n2, uint32_t n3, const std::string s3d_data_file_prefix, const std::string s3d_rdata_file_prefix){
    std::vector<std::string> species = {"H2", "O2", "H2O", "H", "O", "OH"};
    int n_species = species.size();

    std::vector<std::vector<Type>> vars_vec(n_species);
    size_t num_elements = 0;
    for(int i=0; i<n_species; i++){
        auto Xi = MGARD::readfile<Type>((s3d_data_file_prefix + species[i] + ".dat").c_str(), num_elements);
        vars_vec[i] = Xi;
    }

    std::vector<double> value_range(n_species);
    for(int i=0; i<n_species; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }    
    std::vector<double> rel_ebs;
    const int num_snapshot = 18;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_species; i++){
        std::string rdir_prefix = s3d_rdata_file_prefix + species[i];
        for(int j=0; j<num_snapshot; j++){
            std::string filename = rdir_prefix + "_refactored/SZ3_eb_" + std::to_string(j) + ".bin";
            size_t compressed_size = 0;
            auto compressed_data = SZ3_compress_3D(num_elements, n1, n2, n3, vars_vec[i].data(), rel_ebs[j]*value_range[i], compressed_size);
            MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
            free(compressed_data);
        }
    }
}

template<class Type>
void refactor_S3D_SZ3_delta(uint32_t n1, uint32_t n2, uint32_t n3, const std::string s3d_data_file_prefix, const std::string s3d_rdata_file_prefix){
    std::vector<std::string> species = {"H2", "O2", "H2O", "H", "O", "OH"};
    int n_species = species.size();

    std::vector<std::vector<Type>> vars_vec(n_species);
    size_t num_elements = 0;
    for(int i=0; i<n_species; i++){
        auto Xi = MGARD::readfile<Type>((s3d_data_file_prefix + species[i] + ".dat").c_str(), num_elements);
        vars_vec[i] = Xi;
    }

    std::vector<double> value_range(n_species);
    for(int i=0; i<n_species; i++){
        value_range[i] = compute_vr(vars_vec[i]);
        std::cout << "value_range = " << value_range[i] << std::endl;
    }
   
    std::vector<double> rel_ebs;
    const int num_snapshot = 18;
    double eb = 1.0;
    for(int i=0; i<num_snapshot; i++){
        eb /= 10;
        rel_ebs.push_back(eb);
    }
    for(int i=0; i<n_species; i++){
        std::string rdir_prefix = s3d_rdata_file_prefix + species[i];
        std::vector<Type> data_buffer(vars_vec[i]);
        std::vector<Type> dec_data_buffer(vars_vec[i]);
        for(int j=0; j<num_snapshot; j++){
            std::string filename = rdir_prefix + "_refactored/SZ3_delta_eb_" + std::to_string(j) + ".bin";
            size_t compressed_size = 0;
            auto compressed_data = SZ3_compress_3D(num_elements, n1, n2, n3, data_buffer.data(), rel_ebs[j]*value_range[i], compressed_size);
            MGARD::writefile(filename.c_str(), compressed_data, compressed_size);
            SZ3_decompress(compressed_data, compressed_size, dec_data_buffer.data());
            for(int i=0; i<num_elements; i++){
                data_buffer[i] = data_buffer[i] - dec_data_buffer[i];
            }
            free(compressed_data);
        }
    }
}

template<class T>
void refactor_velocities_1D_PMGARD_BP(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    int n_variable = var_list.size();
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());

    uint8_t target_level = 8;   // log2(*min_element(dims.begin(), dims.end())) - 1;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<T>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<T>(decomposer, interleaver, encoder, compressor, collector, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

/* treat 3d data as 1d due to 0-velocity points */
template<class T>
void refactor_velocities_1D_PMGARD_WBP(const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight_for_vtot=4, int max_weight_for_pressure=0, int max_weight_for_density=4, int block_size=1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    std::vector<int> max_weights = {max_weight_for_vtot, max_weight_for_vtot, max_weight_for_vtot, max_weight_for_pressure, max_weight_for_density};
    int n_variable = var_list.size();
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    std::vector<std::vector<T>> Vtot(3, std::vector<T>(num_elements)); // mask not in weight
    std::vector<std::vector<T>> Temp(2, std::vector<T>(num_elements));
    T R = 287.1;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
        if(mask[i]){
            // Vtot[0][i] = fabs(velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[1][i] = fabs(velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[2][i] = fabs(velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // double temp = sqrt(Vtot[0][i]*Vtot[0][i] + Vtot[1][i]*Vtot[1][i] + Vtot[2][i]*Vtot[2][i]);
            // Vtot[0][i] = temp;
            // Vtot[1][i] = temp;
            // Vtot[2][i] = temp;
            // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // double V = velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i];
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        Temp[0][i] = 1.0 / pow(density_vec[i], 2.0);
        Temp[1][i] = 1.0 / pow(density_vec[i], 2.0);
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());

    uint8_t target_level = 8;   // log2(*min_element(dims.begin(), dims.end())) - 1;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto weight_interleaver = MDR::DirectInterleaver<int>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<T>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<T>(decomposer, interleaver, weight_interleaver, encoder, compressor, collector, writer, negabinary);
        if(i < 3) refactor.QoI = Vtot[i];
        else refactor.QoI = Temp[i-3];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weights[i], block_size);            
    }
}

template<class T>
void refactor_velocities_A1D_PMGARD_BP(const std::string data_file_prefix, const std::string rdata_file_prefix, uint32_t n1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());

    uint8_t target_level = 8;   // log2(*min_element(dims.begin(), dims.end())) - 1;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<T>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<T>(decomposer, interleaver, encoder, compressor, collector, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

template<class T>
void refactor_velocities_1D_GE_BP(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    num_valid_data = num_elements;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::GEApproximator<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

template<class T>
void refactor_velocities_1D_GE_WBP(const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight_for_vtot=4, int max_weight_for_pressure=0, int max_weight_for_density=4){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements); 
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    std::vector<int> max_weights = {max_weight_for_vtot, max_weight_for_vtot, max_weight_for_vtot, max_weight_for_pressure, max_weight_for_density};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);

    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    std::vector<std::vector<T>> Vtot(3, std::vector<T>(num_elements));
    std::vector<std::vector<T>> Temp(2, std::vector<T>(num_elements));
    T R = 287.1;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data++;
        }
        if(mask[i]){
            // Vtot[0][i] = fabs(velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[1][i] = fabs(velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[2][i] = fabs(velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // double temp = sqrt(Vtot[0][i]*Vtot[0][i] + Vtot[1][i]*Vtot[1][i] + Vtot[2][i]*Vtot[2][i]);
            // Vtot[0][i] = temp;
            // Vtot[1][i] = temp;
            // Vtot[2][i] = temp;
            // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // double V = velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i];
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        Temp[0][i] = 1.0 / pow(density_vec[i], 2.0);
        Temp[1][i] = 1.0 / pow(density_vec[i], 2.0);
        // if(i == 567082){
        //     std::cout << "index = " << i << ": " << +mask[i] << ", " << Vtot[0][i] << " " << Vtot[1][i] << " " << Vtot[2][i] << std::endl; 
        // }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        std::cout << "metadata location: " << metadata_file << std::endl;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto approximator = PDR::GEApproximator<T>();
        // auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateWBPRefactor_GE<T>(approximator, encoder, compressor, writer, negabinary);
        if(i < 3) refactor.QoI = Vtot[i];
        else refactor.QoI = Temp[i-3];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weights[i]);  
    }
}

template<class T>
void refactor_velocities_1D_Dummy_BP(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    num_valid_data = num_elements;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::DummyApproximator<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

template<class T>
void refactor_velocities_1D_Dummy_WBP(const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight_for_vtot=4, int max_weight_for_pressure=0, int max_weight_for_density=4, int block_size=1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements); 
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    std::vector<int> max_weights = {max_weight_for_vtot, max_weight_for_vtot, max_weight_for_vtot, max_weight_for_pressure, max_weight_for_density};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);

    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    std::vector<std::vector<T>> Vtot(3, std::vector<T>(num_elements));
    std::vector<std::vector<T>> Temp(2, std::vector<T>(num_elements));
    T R = 287.1;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data++;
        }
        if(mask[i]){
            // Vtot[0][i] = fabs(velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[1][i] = fabs(velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[2][i] = fabs(velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // double temp = sqrt(Vtot[0][i]*Vtot[0][i] + Vtot[1][i]*Vtot[1][i] + Vtot[2][i]*Vtot[2][i]);
            // Vtot[0][i] = temp;
            // Vtot[1][i] = temp;
            // Vtot[2][i] = temp;
            // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // double V = velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i];
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        Temp[0][i] = 1.0 / pow(density_vec[i], 2.0);
        Temp[1][i] = 1.0 / pow(density_vec[i], 2.0);
        // if(i == 567082){
        //     std::cout << "index = " << i << ": " << +mask[i] << ", " << Vtot[0][i] << " " << Vtot[1][i] << " " << Vtot[2][i] << std::endl; 
        // }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        std::cout << "metadata location: " << metadata_file << std::endl;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto approximator = PDR::DummyApproximator<T>();
        // auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateWBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        if(i < 3) refactor.QoI = Vtot[i];
        else refactor.QoI = Temp[i-3];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weights[i], block_size);  
    }
}

template<class T>
void refactor_velocities_1D_SZ3_BP(const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    num_valid_data = num_elements;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::SZApproximator<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

template<class T>
void refactor_velocities_1D_SZ3_WBP(const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight_for_vtot=4, int max_weight_for_pressure=0, int max_weight_for_density=4, int block_size=1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    auto pressure_vec = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec, pressure_vec, density_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    std::vector<int> max_weights = {max_weight_for_vtot, max_weight_for_vtot, max_weight_for_vtot, max_weight_for_pressure, max_weight_for_density};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::cout << "num_bitplanes = " << int(num_bitplanes) << std::endl;
    std::vector<uint32_t> dims;
    dims.push_back(num_elements);

    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    std::vector<std::vector<T>> Vtot(3, std::vector<T>(num_elements));
    std::vector<std::vector<T>> Temp(2, std::vector<T>(num_elements));
    T R = 287.1;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data++;
        }
        if(mask[i]){
            // Vtot[0][i] = fabs(velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[1][i] = fabs(velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[2][i] = fabs(velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // double temp = sqrt(Vtot[0][i]*Vtot[0][i] + Vtot[1][i]*Vtot[1][i] + Vtot[2][i]*Vtot[2][i]);
            // Vtot[0][i] = temp;
            // Vtot[1][i] = temp;
            // Vtot[2][i] = temp;
            // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // double V = velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i];
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        Temp[0][i] = 1.0 / pow(density_vec[i], 2.0);
        Temp[1][i] = 1.0 / pow(density_vec[i], 2.0);
        // if(i == 567082){
        //     std::cout << "index = " << i << ": " << +mask[i] << ", " << Vtot[0][i] << " " << Vtot[1][i] << " " << Vtot[2][i] << std::endl; 
        // }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        std::cout << "metadata location: " << metadata_file << std::endl;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateWBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        if(i < 3) refactor.QoI = Vtot[i];
        else refactor.QoI = Temp[i-3];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weights[i], block_size);  
    }
}

template<class T>
void refactor_velocities_A1D_SZ3_BP(const std::string data_file_prefix, const std::string rdata_file_prefix, uint32_t n1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<uint32_t> dims;
    dims.push_back(num_elements);
    // compute masks
    std::vector<unsigned char> mask(num_elements, 0);
    int num_valid_data = 0;
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
            num_valid_data ++;
        }
    }

    std::cout << "num_elements = " << num_elements << ", num_valid_data = " << num_valid_data << std::endl;
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    num_valid_data = num_elements;

    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::SZApproximator<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);            
    }
}

template<class T>
void refactor_velocities_3D_Dummy_BP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims = {n1, n2, n3};

    std::vector<unsigned char> mask(num_elements, 0);
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto approximator = PDR::DummyApproximator<T>();
        // auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);  
    }
}

template<class T>
void refactor_velocities_3D_Dummy_WBP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight=4, int block_size=1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims = {n1, n2, n3};

    std::vector<unsigned char> mask(num_elements, 0);
    std::vector<std::vector<T>> Vtot(n_variable, std::vector<T>(num_elements));
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
        if(mask[i]){
            // Vtot[0][i] = fabs(velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[1][i] = fabs(velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // Vtot[2][i] = fabs(velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]));
            // double temp = sqrt(Vtot[0][i]*Vtot[0][i] + Vtot[1][i]*Vtot[1][i] + Vtot[2][i]*Vtot[2][i]);
            // Vtot[0][i] = temp;
            // Vtot[1][i] = temp;
            // Vtot[2][i] = temp;
            // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // double V = velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i];
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        // if(i == 567082){
        //     std::cout << "index = " << i << ": " << +mask[i] << ", " << Vtot[0][i] << " " << Vtot[1][i] << " " << Vtot[2][i] << std::endl; 
        // }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        std::cout << "metadata location: " << metadata_file << std::endl;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto approximator = PDR::DummyApproximator<T>();
        // auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateWBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.QoI = Vtot[i];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weight, block_size);  
    }
}

template<class T>
void refactor_velocities_3D_SZ3_BP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix){
    std::cout << "refactor_velocities_3D_SZ3_BP" << std::endl;
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims = {n1, n2, n3};

    std::vector<unsigned char> mask(num_elements, 0);
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            std::cout << "filename: " << filename << std::endl;
            files.push_back(filename);
        }
        
        // auto approximator = PDR::DummyApproximator<T>();
        auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);  
    }
}

template<class T>
void refactor_velocities_3D_SZ3_WBP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight=4, int block_size=1){
    std::cout << "refactor_velocities_3D_SZ3_WBP" << std::endl;
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();

    uint8_t target_level = 0;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;
    std::vector<uint32_t> dims = {n1, n2, n3};

    std::vector<unsigned char> mask(num_elements, 0);
    std::vector<std::vector<T>> Vtot(n_variable, std::vector<T>(num_elements));
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
        if(mask[i]){
            // Vtot[0][i] = velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        // auto approximator = PDR::DummyApproximator<T>();
        auto approximator = PDR::SZApproximator<T>();
        // auto encoder = MDR::GroupedBPEncoder<T, uint32_t>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        // auto encoder = MDR::PerBitBPEncoder<T, uint32_t>();
        // negabinary = false;

        // auto compressor = MDR::DefaultLevelCompressor();
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        // auto compressor = MDR::NullLevelCompressor();

        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);
        auto refactor = generateWBPRefactor<T>(approximator, encoder, compressor, writer, negabinary);
        refactor.QoI = Vtot[i];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weight, block_size);  
    }
}

template<class T>
void refactor_velocities_3D_PMGARD_BP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();
    std::vector<uint32_t> dims = {n1, n2, n3};
    uint8_t target_level = 4;// log2(*min_element(dims.begin(), dims.end())) - 1;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<unsigned char> mask(num_elements, 0);
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        
        auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto encoder = MDR::NegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<T>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<T>(decomposer, interleaver, encoder, compressor, collector, writer, negabinary);
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes);  
    }
}

template<class T>
void refactor_velocities_3D_PMGARD_WBP(std::string dataset, uint32_t n1, uint32_t n2, uint32_t n3, const std::string data_file_prefix, const std::string rdata_file_prefix, int max_weight=4, int block_size=1){
    size_t num_elements = 0;
    auto velocityX_vec = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<std::vector<T>> vars_vec = {velocityX_vec, velocityY_vec, velocityZ_vec};
    std::vector<std::string> var_list = {"VelocityX", "VelocityY", "VelocityZ"};
    int n_variable = var_list.size();
    std::vector<uint32_t> dims = {n1, n2, n3};
    uint8_t target_level = 4; // log2(*min_element(dims.begin(), dims.end())) - 1;
    uint8_t num_bitplanes = std::is_same<T, double>::value ? 60 : 32;

    std::vector<unsigned char> mask(num_elements, 0);
    std::vector<std::vector<T>> Vtot(n_variable, std::vector<T>(num_elements));
    for(int i=0; i<num_elements; i++){
        if(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i] != 0){            
            mask[i] = 1;
        }
        if(mask[i]){
            // Vtot[0][i] = velocityX_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[1][i] = velocityY_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            // Vtot[2][i] = velocityZ_vec[i] / std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
            T V = pow(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i], 0.5);
            Vtot[0][i] = 1.0/V;
            Vtot[1][i] = 1.0/V;
            Vtot[2][i] = 1.0/V;
        }
        else{
            Vtot[0][i] = Vtot[1][i] = Vtot[2][i] = 0;
        }
        // Vtot[0][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[1][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
        // Vtot[2][i] = std::sqrt(velocityX_vec[i]*velocityX_vec[i] + velocityY_vec[i]*velocityY_vec[i] + velocityZ_vec[i]*velocityZ_vec[i]);
    }
    std::string mask_file = rdata_file_prefix + "mask.bin";
    MGARD::writefile(mask_file.c_str(), mask.data(), mask.size());
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_list[i];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = target_level + 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto weight_interleaver = MDR::DirectInterleaver<int>();
        auto encoder = MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>();
        bool negabinary = true;
        auto compressor = MDR::AdaptiveLevelCompressor(64);
        auto collector = MDR::SquaredErrorCollector<T>();
        auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
        auto refactor = generateRefactor<T>(decomposer, interleaver, weight_interleaver, encoder, compressor, collector, writer, negabinary);
        refactor.QoI = Vtot[i];
        refactor.refactor(vars_vec[i].data(), dims, target_level, num_bitplanes, max_weight, block_size);  
    }
}
}
#endif