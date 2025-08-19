#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"
#include "nomask_Synthesizer4GE.hpp"
#include "PDR/Reconstructor/Reconstructor.hpp"
#define Dummy 0
#define SZ3 1
#define PMGARD 2
#define GE 3
using namespace MDR;

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
double * mu_ori = NULL;
std::vector<double> error_mu;
std::vector<double> error_est_mu;
int iter = 0;

template<class T>
bool halfing_error_mu_uniform(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	T c_3 = T_r + S;
	T max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		T TrS_TS = c_3 / (Temp + S);
		T e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		T T_Tr_3 = pow(Temp/T_r, 3);
		T e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		T T_Tr_3_sqrt = sqrt(T_Tr_3);
		T e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		T mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_P = ebs[0];
		T eb_D = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			T Temp = P[i] / (D[i] * R);
			T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
			double TrS_TS = c_3 / (Temp + S);
			double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
			double T_Tr_3 = pow(Temp/T_r, 3);
			double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
			double T_Tr_3_sqrt = sqrt(T_Tr_3);
			estimate_error = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
            if (ebs[0] / eb_P > 10) break;
		}
		ebs[0] = eb_P;
		ebs[1] = eb_D;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_mu_coordinate(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	T c_3 = T_r + S;
	T max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		T Temp = P[i] / (D[i] * R);
		T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		T TrS_TS = c_3 / (Temp + S);
		T e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		T T_Tr_3 = pow(Temp/T_r, 3);
		T e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		T T_Tr_3_sqrt = sqrt(T_Tr_3);
		T e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		T mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;

    // ************************************************ P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_P = ebs[0];
		T eb_D = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            T estimate_error_P = 0;
			{
                T eb_P_ = eb_P / 1.5;
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P_, eb_D);
                T Temp = P[i] / (D[i] * R);
                T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
                double TrS_TS = c_3 / (Temp + S);
                double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
                double T_Tr_3 = pow(Temp/T_r, 3);
                double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
                double T_Tr_3_sqrt = sqrt(T_Tr_3);
                estimate_error_P = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
            }
            T estimate_error_D = 0;
			{
                T eb_D_ = eb_D / 1.5;
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D_);
                T Temp = P[i] / (D[i] * R);
                T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
                double TrS_TS = c_3 / (Temp + S);
                double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
                double T_Tr_3 = pow(Temp/T_r, 3);
                double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
                double T_Tr_3_sqrt = sqrt(T_Tr_3);
                estimate_error_D = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
            }
            std::cout << estimate_error_P << " " << estimate_error_D << std::endl;
            const T epsilon = 1e-6;
            T min_error = std::min({estimate_error_P, estimate_error_D});
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
bool halfing_error_mu_uniform(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	T c_3 = T_r + S;
	T max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
		T Temp = P[i] / (D[i] * R);
		T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		T TrS_TS = c_3 / (Temp + S);
		T e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		T T_Tr_3 = pow(Temp/T_r, 3);
		T e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		T T_Tr_3_sqrt = sqrt(T_Tr_3);
		T e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		T mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_P = ebs[0];
		T eb_D = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
			T Temp = P[i] / (D[i] * R);
			T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
			double TrS_TS = c_3 / (Temp + S);
			double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
			double T_Tr_3 = pow(Temp/T_r, 3);
			double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
			double T_Tr_3_sqrt = sqrt(T_Tr_3);
			estimate_error = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
            if (ebs[0] / eb_P > 10) break;
		}
		ebs[0] = eb_P;
		ebs[1] = eb_D;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_mu_coordinate(const T * P, const T * D, size_t n, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_P = ebs[0];
	T eb_D = ebs[1];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	T c_3 = T_r + S;
	T max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
		T Temp = P[i] / (D[i] * R);
		T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		T TrS_TS = c_3 / (Temp + S);
		T e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		T T_Tr_3 = pow(Temp/T_r, 3);
		T e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		T T_Tr_3_sqrt = sqrt(T_Tr_3);
		T e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		T mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;

    // ************************************************ P and D ************************************************
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_P = ebs[0];
		T eb_D = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            T estimate_error_P = 0;
			{
                T eb_P_ = eb_P / 1.5;
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P_ / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
                T Temp = P[i] / (D[i] * R);
                T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
                double TrS_TS = c_3 / (Temp + S);
                double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
                double T_Tr_3 = pow(Temp/T_r, 3);
                double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
                double T_Tr_3_sqrt = sqrt(T_Tr_3);
                estimate_error_P = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
            }
            T estimate_error_D = 0;
			{
                T eb_D_ = eb_D / 1.5;
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D_ / static_cast<T>(std::pow(2.0, weights[1][i])));
                T Temp = P[i] / (D[i] * R);
                T e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
                double TrS_TS = c_3 / (Temp + S);
                double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
                double T_Tr_3 = pow(Temp/T_r, 3);
                double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
                double T_Tr_3_sqrt = sqrt(T_Tr_3);
                estimate_error_D = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
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
std::vector<size_t> retrieve_mu_Dummy(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_mu.data()), error_est_mu.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_mu_SZ3(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
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
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs);
            else tolerance_met = halfing_error_mu_coordinate(P_dec, D_dec, num_elements, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            else tolerance_met = halfing_error_mu_coordinate(P_dec, D_dec, num_elements, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_mu.data()), error_est_mu.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_mu_PMGARD(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
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
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
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
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            MDR::print_vec(ebs);
            tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_mu.data()), error_est_mu.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_mu_GE(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::GEApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = NegaBinaryBPEncoder<T, uint32_t>();
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
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs);
            else tolerance_met = halfing_error_mu_coordinate(P_dec, D_dec, num_elements, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    else{
        std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
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
            error_mu = std::vector<T>(num_elements);
            error_est_mu = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
            else tolerance_met = halfing_error_mu_coordinate(P_dec, D_dec, num_elements, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_mu.data()), error_est_mu.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_mu);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_mu);   	
        }
    }
    return total_retrieved_size;
}

int main(int argc, char ** argv){

    using T = double;
	int argv_id = 1;
    int compressor = atoi(argv[argv_id++]);
    int weighted = atoi(argv[argv_id++]);
    int decrease_method = atoi(argv[argv_id++]);
    T target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<T> ebs;
    ebs.push_back(compute_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(D_ori)*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(P_ori.data(), P_ori.size())*target_rel_eb);
    // ebs.push_back(compute_max_abs_value(D_ori.data(), D_ori.size())*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {P_ori, D_ori};
    std::vector<double> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> mu(num_elements);
    compute_mu(P_ori.data(), D_ori.data(), num_elements, mu.data());
	mu_ori = mu.data();
    T tau = compute_value_range(mu)*target_rel_eb;
    T max_act_error = 0, max_est_error = 0;
    size_t weight_file_size = 0;
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    switch (compressor)
    {
    case Dummy:
        total_retrieved_size = retrieve_mu_Dummy<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
        break;
    case SZ3:
        total_retrieved_size = retrieve_mu_SZ3<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
        break;
    case PMGARD:
        total_retrieved_size = retrieve_mu_PMGARD<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
        break;
    case GE:
        total_retrieved_size = retrieve_mu_GE<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
        break;
    default:
        break;
    }

	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested_error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
	std::cout << "iter = " << iter << std::endl;
   
   	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), size_t(0)) + weight_file_size;
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_size[i] << ", ";
    }
    std::cout << "weight_file_size = " << weight_file_size << std::endl;
	// MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;
    std::cout << "bitrate = " << ((sizeof(T) * 8) / cr) << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}