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
#define Dummy_Cmp 0
#define SZ3_Cmp 1
#define PMGARD 2
#define GE_Cmp 3
#define HPEZ_Cmp 4
#define SZ2_Cmp 5
#define MGARD_Cmp 6
#define HPEZ_Vtot2 7
using namespace MDR;

using T = double;
using T_stream = uint32_t;
std::vector<double> Vx_ori;
std::vector<double> Vy_ori;
std::vector<double> Vz_ori;
double * Vx_dec = NULL;
double * Vy_dec = NULL;
double * Vz_dec = NULL;
double * V_TOT2_ori = NULL;
std::vector<double> error_V_TOT2;
std::vector<double> error_est_V_TOT2;
int iter = 0;


template<class T>
bool halfing_error_V_TOT2_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T max_value = 0;
	int max_index = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		// T e_V_TOT = 0;
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT2[i] = e_V_TOT_2;
		error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

		if(max_value < error_est_V_TOT2[i]){
			max_value = error_est_V_TOT2[i];
			max_index = i;
			// max_weight_index = weight_index;
		}
		// if(mask[i]) weight_index++;
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		// estimate
		auto i = max_index;
		T estimate_error = max_value;
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
    		// std::cout << "uniform decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5;		        		
			estimate_error = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT2_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T max_value = 0;
	int max_index = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT2[i] = e_V_TOT_2;
		error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

		if(max_value < error_est_V_TOT2[i]){
			max_value = error_est_V_TOT2[i];
			max_index = i;
			// max_weight_index = weight_index;
		}
		// if(mask[i]) weight_index++;
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;

	// ************************************************ Vx, Vy and Vz ************************************************
	if(max_value > tau){
		// estimate
		auto i = max_index;
		T estimate_error = max_value;
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		while(estimate_error > tau){
			// std::cout << "coordinate decrease\n";
			// std::cout << "coordinate decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
    		T estimate_error_Vx = 0;
			{
				T eb_Vx_ = eb_Vx / 1.5;
				estimate_error_Vx = compute_bound_x_square(Vx[i], eb_Vx_) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
			}
			T estimate_error_Vy = 0;
			{
				T eb_Vy_ = eb_Vy / 1.5;
				estimate_error_Vy = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy_) + compute_bound_x_square(Vz[i], eb_Vz);
			}
			T estimate_error_Vz = 0;
			{
				T eb_Vz_ = eb_Vz / 1.5;
				estimate_error_Vz = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz_);
			}
			// // std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << std::endl;
			// const T epsilon = 1e-6;
            // T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz});
            // bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            // bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            // bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            // estimate_error = min_error;
            // if (close_Vx) eb_Vx /= 1.5;
            // if (close_Vy) eb_Vy /= 1.5;
            // if (close_Vz) eb_Vz /= 1.5;
			// if (ebs[0] / eb_Vx > 10 || ebs[1] / eb_Vy > 10 || ebs[2] / eb_Vz > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT2_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T max_value = 0;
	int max_index = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		if(mask[i]) {
			// error of total velocity square
			T e_V_TOT_2 = 0;
			e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
			T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

			error_est_V_TOT2[i] = e_V_TOT_2;
			error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

			if(max_value < error_est_V_TOT2[i]){
				max_value = error_est_V_TOT2[i];
				max_index = i;
				// max_weight_index = weight_index;
			}
			// if(mask[i]) weight_index++;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		// estimate
		auto i = max_index;
		T estimate_error = max_value;
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5;		        		
			estimate_error = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT2_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T max_value = 0;
	int max_index = 0;
	for(int i=0; i<n; i++){
		if(mask[i]) {
			// error of total velocity square
			T e_V_TOT_2 = 0;
			e_V_TOT_2 = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
			T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

			error_est_V_TOT2[i] = e_V_TOT_2;
			error_V_TOT2[i] = V_TOT_2 - V_TOT2_ori[i];

			if(max_value < error_est_V_TOT2[i]){
				max_value = error_est_V_TOT2[i];
				max_index = i;
				// max_weight_index = weight_index;
			}
			// if(mask[i]) weight_index++;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;

	// ************************************************ Vx, Vy and Vz ************************************************
	if(max_value > tau){
		// estimate
		auto i = max_index;
		T estimate_error = max_value;
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		while(estimate_error > tau){
			// std::cout << "coordinate decrease\n";
			// std::cout << "coordinate decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
    		T estimate_error_Vx = 0;
			{
				T eb_Vx_ = eb_Vx / 1.5;
				estimate_error_Vx = compute_bound_x_square(Vx[i], ldexp(eb_Vx_, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
			}
			T estimate_error_Vy = 0;
			{
				T eb_Vy_ = eb_Vy / 1.5;
				estimate_error_Vy = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy_, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
			}
			T estimate_error_Vz = 0;
			{
				T eb_Vz_ = eb_Vz / 1.5;
				estimate_error_Vz = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz_, -weights[2][i]));
			}
			// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << std::endl;
			T sum_err = 3 * estimate_error - (estimate_error_Vx + estimate_error_Vy + estimate_error_Vz);
            T w_Vx = (estimate_error - estimate_error_Vx) / sum_err;
            T w_Vy = (estimate_error - estimate_error_Vy) / sum_err;
            T w_Vz = (estimate_error - estimate_error_Vz) / sum_err;
            // === Smooth proportional update ===
            T factor_base = 1.5;
            T alpha = 1.0 - (1.0 / factor_base);
            eb_Vx *= (1.0 - alpha * w_Vx);
            eb_Vy *= (1.0 - alpha * w_Vy);
            eb_Vz *= (1.0 - alpha * w_Vz);

			{
				estimate_error = compute_bound_x_square(Vx[i], ldexp(eb_Vx, -weights[0][i])) + compute_bound_x_square(Vy[i], ldexp(eb_Vy, -weights[1][i])) + compute_bound_x_square(Vz[i], ldexp(eb_Vz, -weights[2][i]));
			}
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_Dummy(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
    	}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(T));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::DummyApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
			reconstructors.back().load_weight();
			reconstructors.back().span_weight();
			weights[i] = reconstructors.back().get_int_weights();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_SZ3(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZ3Approximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::SZ3Approximator<T>();
			auto encoder = NegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(T));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZ3Approximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::SZ3Approximator<T>();
			auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
			weight_file_size += reconstructors.back().get_weight_file_size();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_PMGARD(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
	int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				if(i < 3){
                    // reconstruct with mask
                    int index = 0;
                    for(int j=0; j<num_elements; j++){
                        if(mask[j]){
                            reconstructed_vars[i][j] = reconstructed_data[index ++];
                        }
                        else{
                            reconstructed_vars[i][j] = 0; 
                            index++;
                        }
                    }
                }
                else{
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<MDR::WeightReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, DirectInterleaver<int>, WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
				// auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				if(i < 3){
                    // reconstruct with mask
                    int index = 0;
                    for(int j=0; j<num_elements; j++){
                        if(mask[j]){
                            reconstructed_vars[i][j] = reconstructed_data[index ++];
                        }
                        else{
                            reconstructed_vars[i][j] = 0; 
                            index++;
                        }
                    }
                }
                else{
                    memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                }
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_GE(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::GEApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
    	}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(T));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
			reconstructors.back().mask = mask;
			reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
			weight_file_size += reconstructors.back().get_weight_file_size();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_HPEZ(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::HPEZApproximator<T>();
			auto encoder = NegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			// std::cout << "iter " << iter << std::endl;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(T));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::HPEZApproximator<T>();
			auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			if(i==0) reconstructors.back().fetch_weight = true;
			else reconstructors.back().copy_int_weights(weights[0]);
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[0].get_max_weight());
		}
		weight_file_size = reconstructors[0].get_weight_file_size();
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			// std::cout << "iter " << iter << std::endl;
			for(int i=0; i<n_variable; i++){
				// Timer reconstructor_timer;
				// reconstructor_timer.start();
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[0].get_max_weight()), -1);
				// reconstructor_timer.end();
				// reconstructor_timer.print("Reconstructor");
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			// Timer decrease_timer;
			// decrease_timer.start();
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// decrease_timer.end();
			// decrease_timer.print("Coordinate Decrease");
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_SZ2(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZ2Approximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::SZ2Approximator<T>();
			auto encoder = NegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZ2Approximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::SZ2Approximator<T>();
			auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
			weight_file_size += reconstructors.back().get_weight_file_size();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_MGARD(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::MGARDApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::MGARDApproximator<T>();
			auto encoder = NegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::MGARDApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::MGARDApproximator<T>();
			auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
			weight_file_size += reconstructors.back().get_weight_file_size();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT2_HPEZ_square(std::string rdata_file_prefix, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::HPEZApproximator<T>();
			auto encoder = NegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().load_metadata();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			// std::cout << "iter " << iter << std::endl;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(T));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements));
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
			std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
			std::vector<std::string> files;
			int num_levels = 1;
			for(int i=0; i<num_levels; i++){
				std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
				files.push_back(filename);
			}
			auto approximator = PDR::HPEZApproximator<T>();
			auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
			auto compressor = AdaptiveLevelCompressor(64);
			auto estimator = MaxErrorEstimatorHB<T>();
			auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
			auto retriever = ConcatLevelFileRetriever(metadata_file, files);
			reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
			reconstructors.back().mask = mask;
			reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
			ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
			weight_file_size += reconstructors.back().get_weight_file_size();
		}
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			// std::cout << "iter " << iter << std::endl;
			for(int i=0; i<n_variable; i++){
				// Timer reconstructor_timer;
				// reconstructor_timer.start();
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
				// reconstructor_timer.end();
				// reconstructor_timer.print("Reconstructor");
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			// for(int i=0; i<n_variable; i++){
			// 	std::cout << total_retrieved_size[i] << ", ";
			// }
			// std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			// MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			// MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			// MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT2 = std::vector<T>(num_elements);
			error_est_V_TOT2 = std::vector<T>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			// Timer decrease_timer;
			// decrease_timer.start();
			if(!decrease_method) tolerance_met = halfing_error_V_TOT2_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT2_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs, weights);
			// decrease_timer.end();
			// decrease_timer.print("Coordinate Decrease");
			// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			/* test
			std::string filename = "./Result/Vtot_err.dat";
			std::ofstream outfile1(filename, std::ios::binary);
			if (!outfile1.is_open()) {
				std::cerr << "Failed to open file for writing: " << filename << std::endl;
				exit(-1);
			}

			outfile1.write(reinterpret_cast<const char*>(error_est_V_TOT.data()), error_est_V_TOT.size() * sizeof(float));
		
			outfile1.close();
			std::cout << "Data saved successfully to " << filename << std::endl;
			//*/
			// std::cout << names[0] << " requested error = " << tau << std::endl;
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT2);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT2);  
		}
	}
	return total_retrieved_size;
}

int main(int argc, char ** argv){
	using T = double;
	using T_stream = uint32_t;
	int argv_id = 1;
	int compressor = atoi(argv[argv_id++]);
    int weighted = atoi(argv[argv_id++]);
	int decrease_method = atoi(argv[argv_id++]);
    T target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<T> ebs;
	ebs.push_back(compute_value_range(Vx_ori)*target_rel_eb);
	ebs.push_back(compute_value_range(Vy_ori)*target_rel_eb);
	ebs.push_back(compute_value_range(Vz_ori)*target_rel_eb);
	// ebs.push_back(compute_max_abs_value(Vx_ori.data(), Vx_ori.size()*target_rel_eb));
	// ebs.push_back(compute_max_abs_value(Vy_ori.data(), Vy_ori.size()*target_rel_eb));
	// ebs.push_back(compute_max_abs_value(Vz_ori.data(), Vz_ori.size()*target_rel_eb));
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Vx_ori, Vy_ori, Vz_ori};
    std::vector<T> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> V_TOT2(num_elements);
    compute_VTOT2(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT2.data());
	V_TOT2_ori = V_TOT2.data();
	/* test
    std::string filename = "./Result/Vtot_ori.dat";
    std::ofstream outfile1(filename, std::ios::binary);
    if (!outfile1.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
		exit(-1);
    }
	std::cout << "VTOT.size: " << V_TOT.size() << std::endl;
    outfile1.write(reinterpret_cast<const char*>(V_TOT.data()), V_TOT.size() * sizeof(T));
    
    outfile1.close();
    std::cout << "Data saved successfully to " << filename << std::endl;
    //*/
	// std::cout << "Value Range of V_TOT = " << compute_value_range(V_TOT) << std::endl;
    T tau = compute_value_range(V_TOT2)*target_rel_eb;
	// std::cout << "Tau = " << tau << std::endl;

    std::string mask_file = rdata_file_prefix + "mask.bin";
    uint32_t mask_file_size = 0;
    auto mask = readmask(mask_file.c_str(), mask_file_size);
	T max_act_error = 0;
	T max_est_error = 0;
	size_t weight_file_size = 0;
	std::vector<size_t> total_retrieved_size;
	switch (compressor)
	{
	case Dummy_Cmp:
		total_retrieved_size = retrieve_V_TOT2_Dummy<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case SZ3_Cmp:
		total_retrieved_size = retrieve_V_TOT2_SZ3<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case PMGARD:
		total_retrieved_size = retrieve_V_TOT2_PMGARD<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case GE_Cmp:
		total_retrieved_size = retrieve_V_TOT2_GE<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case HPEZ_Cmp:
		total_retrieved_size = retrieve_V_TOT2_HPEZ<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case SZ2_Cmp:
		total_retrieved_size = retrieve_V_TOT2_SZ2<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case MGARD_Cmp:
		total_retrieved_size = retrieve_V_TOT2_MGARD<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
		break;
	case HPEZ_Vtot2:
		total_retrieved_size = retrieve_V_TOT2_HPEZ_square<T>(rdata_file_prefix, tau, ebs, num_elements, mask, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
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

	size_t total_size = mask_file_size + std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), size_t(0)) + weight_file_size;
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_size[i] << ", ";
    }
	std::cout << "mask_file_size = " << mask_file_size << ", weight_file_size = " << weight_file_size << std::endl;
	// MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;
	std::cout << "bitrate = " << ((sizeof(T) * 8) / cr) << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}