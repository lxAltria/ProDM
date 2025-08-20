#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include<limits>
#include "mpi.h"
#include <sstream>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "nomask_Synthesizer4GE.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#include "MDR/RefactorUtils.hpp"
#define Dummy 0
#define SZ3 1
#define PMGARD 2
#define GE 3
using namespace MDR;

using T = double;
using T_stream = uint32_t;
std::vector<double> P_ori;
std::vector<double> D_ori;
std::vector<double> Vx_ori;
std::vector<double> Vy_ori;
std::vector<double> Vz_ori;
double * P_dec = NULL;
double * D_dec = NULL;
double * Vx_dec = NULL;
double * Vy_dec = NULL;
double * Vz_dec = NULL;
double * Temp_ori = NULL;
std::vector<double> error_Temp;
std::vector<double> error_est_Temp;
int iter = 0;
double local_elapsed_time = 0;

template<class T>
bool halfing_error_T_uniform(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T c_1 = 1.0 / R;
	T max_value = 0;;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);

		error_est_Temp[i] = e_T;
		error_Temp[i] = Temp - Temp_ori[i];

		if(max_value < error_est_Temp[i]){
			max_value = error_est_Temp[i];
			max_index = i;
		}
	}
	// std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	// std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	// std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	// std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
        // float T = c_1 * P[i] / D[i];
        T eb_P = ebs[0];
        T eb_D = ebs[1];
        while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
            eb_P = eb_P / 1.5;
            eb_D = eb_D / 1.5;
            estimate_error = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
            if (ebs[0] / eb_P > 10) break;
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_T_coordinate(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T c_1 = 1.0 / R;
	T max_value = 0;;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);

		error_est_Temp[i] = e_T;
		error_Temp[i] = Temp - Temp_ori[i];

		if(max_value < error_est_Temp[i]){
			max_value = error_est_Temp[i];
			max_index = i;
		}
	}
	// std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	// std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	// std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	// std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;

    // ************************************************ P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
        // float T = c_1 * P[i] / D[i];
        T eb_P = ebs[0];
        T eb_D = ebs[1];
        while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            T estimate_error_P = 0;
            {
                T eb_P_ = eb_P / 1.5;
                estimate_error_P = c_1 * compute_bound_division(P[i], D[i], eb_P_, eb_D);
            }
            T estimate_error_D = 0;
            {
                T eb_D_ = eb_D / 1.5;
                estimate_error_D = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D_);
            }
            // std::cout << estimate_error_P << " " << estimate_error_D << std::endl;
            const T relative_epsilon = 1e-3;
            T min_error = std::min({estimate_error_P, estimate_error_D});
            T epsilon = std::max(relative_epsilon * min_error, static_cast<T>(1e-12));
            bool close_P  = fabs(estimate_error_P - min_error) < epsilon;
            bool close_D  = fabs(estimate_error_D - min_error) < epsilon;
            estimate_error = min_error;
            if (close_P)  eb_P /= 1.5;
            if (close_D)  eb_D /= 1.5;
            if (ebs[0] / eb_P > 10 || ebs[1] / eb_D > 10) break;
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}

	return true;
}

template<class T>
bool halfing_error_T_uniform(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T c_1 = 1.0 / R;
	T max_value = 0;;
	int max_index = 0;
    T max_e_T = 0;
	T max_T = 0;
	T max_P = 0;
	T max_D = 0;
	int max_weight_P = 0;
	int max_weight_D = 0;
	for(int i=0; i<n; i++){
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);

		error_est_Temp[i] = e_T;
		error_Temp[i] = Temp - Temp_ori[i];

		if(max_value < error_est_Temp[i]){
			max_value = error_est_Temp[i];
			max_index = i;
            max_e_T = e_T;
            max_T = Temp;
            max_P = P[i];
            max_D = D[i];
            max_weight_P = weights[0][i];
            max_weight_D = weights[1][i];
		}
	}
	// std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	// std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	// std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	// std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    // std::cout << "max_weight_P = " << max_weight_P << ", max_weight_D = " << max_weight_D << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
        // float T = c_1 * P[i] / D[i];
        T eb_P = ebs[0];
        T eb_D = ebs[1];
        while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
            eb_P = eb_P / 1.5;
            eb_D = eb_D / 1.5;
            estimate_error = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
            if (ebs[0] / eb_P > 10) break;
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}
	return true;
} 

template<class T>
bool halfing_error_T_coordinate(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T c_1 = 1.0 / R;
	T max_value = 0;;
	int max_index = 0;
    T max_e_T = 0;
	T max_T = 0;
	T max_P = 0;
	T max_D = 0;
	int max_weight_P = 0;
	int max_weight_D = 0;
	for(int i=0; i<n; i++){
		// error of temperature
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
		T Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);

		error_est_Temp[i] = e_T;
		error_Temp[i] = Temp - Temp_ori[i];

		if(max_value < error_est_Temp[i]){
			max_value = error_est_Temp[i];
			max_index = i;
            max_e_T = e_T;
            max_T = Temp;
            max_P = P[i];
            max_D = D[i];
            max_weight_P = weights[0][i];
            max_weight_D = weights[1][i];
		}
	}
	// std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	// std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	// std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	// std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    // std::cout << "max_weight_P = " << max_weight_P << ", max_weight_D = " << max_weight_D << std::endl;

    // ************************************************ P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
        // float T = c_1 * P[i] / D[i];
        T eb_P = ebs[0];
        T eb_D = ebs[1];
        while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            T estimate_error_P = 0;
            {
                T eb_P_ = eb_P / 1.5;
                estimate_error_P = c_1 * compute_bound_division(P[i], D[i], eb_P_ / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
            }
            T estimate_error_D = 0;
            {
                T eb_D_ = eb_D / 1.5;
                estimate_error_D = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D_ / static_cast<T>(std::pow(2.0, weights[1][i])));
            }
            // std::cout << estimate_error_P << " " << estimate_error_D << std::endl;
            const T relative_epsilon = 1e-3;
            T min_error = std::min({estimate_error_P, estimate_error_D});
            T epsilon = std::max(relative_epsilon * min_error, static_cast<T>(1e-12));
            bool close_P  = fabs(estimate_error_P - min_error) < epsilon;
            bool close_D  = fabs(estimate_error_D - min_error) < epsilon;
            estimate_error = min_error;
            if (close_P)  eb_P /= 1.5;
            if (close_D)  eb_D /= 1.5;
            if (ebs[0] / eb_P > 10 || ebs[1] / eb_D > 10) break;
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}

	return true;
} 

template<class T>
T compute_global_value_range(const std::vector<T> data_vec){
	T global_max = 0, global_min = 0;
	T local_max = -std::numeric_limits<T>::max();
	T local_min = std::numeric_limits<T>::max();
	for(int i=0; i<data_vec.size(); i++){
		if(data_vec[i] > local_max) local_max = data_vec[i];
		if(data_vec[i] < local_min)	local_min = data_vec[i];
	}
	if(std::is_same<T, double>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	else if(std::is_same<T, float>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
	}
	return global_max - global_min;
}

template<class T>
std::vector<size_t> retrieve_T_Dummy(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::DummyApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
        }
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::DummyApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
            reconstructors.back().load_weight();
            reconstructors.back().span_weight();
            weights[i] = reconstructors.back().get_int_weights();
        }
        weight_file_size = reconstructors[n_variable - 1].get_weight_file_size();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Temp.data()), error_est_Temp.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_T_SZ3(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::SZApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
        }
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
            else tolerance_met = halfing_error_T_coordinate(P_dec, D_dec, num_elements, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::SZApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i==0) reconstructors.back().fetch_weight = true;
			else reconstructors.back().copy_int_weights(weights[0]);
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
            ebs[i] *= static_cast<T>(std::pow(2.0, reconstructors[0].get_max_weight()));
        }    
        weight_file_size = reconstructors[0].get_weight_file_size();
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[0].get_max_weight())), -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            else tolerance_met = halfing_error_T_coordinate(P_dec, D_dec, num_elements, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Temp.data()), error_est_Temp.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_T_PMGARD(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 9;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = DirectInterleaver<T>();
            auto encoder = NegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
        }
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
    }
    else{
        std::vector<MDR::WeightReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, DirectInterleaver<int>, WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 9;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = DirectInterleaver<T>();
            auto weight_interleaver = DirectInterleaver<int>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, weight_interleaver, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
            reconstructors.back().load_weight();
            reconstructors.back().span_weight();
            weights[i] = reconstructors.back().get_int_weights();
        }
        weight_file_size = reconstructors[n_variable - 1].get_weight_file_size();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Temp.data()), error_est_Temp.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_T_GE(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::GEApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::GEApproximator<T>();
            auto encoder = NegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().load_metadata();
        }
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
            else tolerance_met = halfing_error_T_coordinate(P_dec, D_dec, num_elements, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    else{
        std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::GEApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor_GE<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            if(i==0) reconstructors.back().fetch_weight = true;
			else reconstructors.back().copy_int_weights(weights[0]);
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
            ebs[i] *= static_cast<T>(std::pow(2.0, reconstructors[0].get_max_weight()));
        }
        weight_file_size = reconstructors[0].get_weight_file_size();
        local_elapsed_time = -MPI_Wtime();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[0].get_max_weight())), -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            P_dec = reconstructed_vars[0].data();
            D_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_Temp = std::vector<T>(num_elements);
            error_est_Temp = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            else tolerance_met = halfing_error_T_coordinate(P_dec, D_dec, num_elements, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_Temp.data()), error_est_Temp.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_Temp);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
        }
        local_elapsed_time += MPI_Wtime();
    }
    return total_retrieved_size;
}

int main(int argc, char ** argv){
    if(argc != 7){
		std::cout << "usage: para_VTOT [Approximator] [weighted] [decrease_method] [eb] [input path] [output path]" << std::endl;
 	}
	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	std::ostringstream oss;
	oss << rank;

    using T = double;
	int argv_id = 1;
    int compressor = atoi(argv[argv_id++]);
    int weighted = atoi(argv[argv_id++]);
    int decrease_method = atoi(argv[argv_id++]);
    double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
    std::string output_path = argv[argv_id++];

    data_prefix_path += oss.str();
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<T> ebs;
    ebs.push_back(compute_global_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(D_ori)*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(P_ori.data(), P_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(D_ori.data(), D_ori.size())*target_rel_eb);
	int n_variable = ebs.size();

    std::vector<T> Temp(num_elements);
    compute_T(P_ori.data(), D_ori.data(), num_elements, Temp.data());

    double tau, global_tau;
    // if(rank == 0) std::cout << "ebs[0] = " << ebs[0] << ", ebs[1] = " << ebs[1] << std::endl;
	global_tau = compute_global_value_range(Temp) * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
	Temp_ori = Temp.data();
    
    T local_max_act_error = 0;
	T local_max_est_error = 0;
    size_t weight_file_size = 0;
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    switch (compressor)
    {
    case Dummy:
        total_retrieved_size = retrieve_T_Dummy<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
        break;
    case SZ3:
        total_retrieved_size = retrieve_T_SZ3<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
        break;
    case PMGARD:
        total_retrieved_size = retrieve_T_PMGARD<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
        break;
    case GE:
        total_retrieved_size = retrieve_T_GE<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
        break;
    default:
        break;
    }
    double global_elapsed_time = 0;
    MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;

    if(!rank) printf("requested_error = %.10f\n", global_tau);
	
    double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_est_error = %.10f\n", global_max_est_error);

	double global_max_act_error = 0;
	MPI_Reduce(&local_max_act_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_act_error = %.10f\n", global_max_act_error);

	unsigned long long local_total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0ULL) + static_cast<unsigned long long>(weight_file_size);

    // for(int i=0; i<size; i++){
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	if(rank == i) std::cout << "Rank " << rank << ": local_bitrate = " << 8 * local_total_size * 1.0 / (num_elements * n_variable) << std::endl;
	// }

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);

    MPI_Finalize();

    return 0;
}