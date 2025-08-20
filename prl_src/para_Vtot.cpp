#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <limits>
#include "mpi.h"
#include <sstream>
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

const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ"};

using T = float;
using T_stream = uint32_t;
std::vector<float> P_ori;
std::vector<float> D_ori;
std::vector<float> Vx_ori;
std::vector<float> Vy_ori;
std::vector<float> Vz_ori;
float * P_dec = NULL;
float * D_dec = NULL;
float * Vx_dec = NULL;
float * Vy_dec = NULL;
float * Vz_dec = NULL;
float * V_TOT_ori = NULL;
std::vector<double> error_V_TOT;
std::vector<double> error_est_V_TOT;
int iter = 0;
double local_elapsed_time = 0;

template<class T>
bool halfing_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	double max_e_V_TOT_2 = 0;
	double max_V_TOT_2 = 0;
	double max_Vx = 0;
	double max_Vy = 0;
	double max_Vz = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_e_V_TOT_2 = e_V_TOT_2;
			max_V_TOT_2 = V_TOT_2;
			max_Vx = Vx[i];
			max_Vy = Vy[i];
			max_Vz = Vz[i];
			max_index = i;
			// max_weight_index = weight_index;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		// estimate
		auto i = max_index;
		double estimate_error = max_value;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
    		// std::cout << "uniform decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5;		        		
			double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
			// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			if (ebs[0] / eb_Vx > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	double max_e_V_TOT_2 = 0;
	double max_V_TOT_2 = 0;
	double max_Vx = 0;
	double max_Vy = 0;
	double max_Vz = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_e_V_TOT_2 = e_V_TOT_2;
			max_V_TOT_2 = V_TOT_2;
			max_Vx = Vx[i];
			max_Vy = Vy[i];
			max_Vz = Vz[i];
			max_index = i;
			// max_weight_index = weight_index;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;

	// ************************************************ Vx, Vy and Vz ************************************************
	if(max_value > tau){
		// estimate
		auto i = max_index;
		double estimate_error = max_value;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		while(estimate_error > tau){
			// std::cout << "coordinate decrease\n";
			// std::cout << "coordinate decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
    		double estimate_error_Vx = 0;
			{
				double eb_Vx_ = eb_Vx / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx_) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vx = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			double estimate_error_Vy = 0;
			{
				double eb_Vy_ = eb_Vy / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy_) + compute_bound_x_square((double) Vz[i], eb_Vz);
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vy = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			double estimate_error_Vz = 0;
			{
				double eb_Vz_ = eb_Vz / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz_);
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vz = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << std::endl;
			const double relative_epsilon = 1e-3;
			double min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz});
			double epsilon = std::max(relative_epsilon * min_error, static_cast<double>(1e-12));
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
			if (ebs[0] / eb_Vx > 10 || ebs[1] / eb_Vy > 10 || ebs[2] / eb_Vz > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const double tau, std::vector<double>& ebs, std::vector<std::vector<int>> weights){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	double max_e_V_TOT_2 = 0;
	double max_V_TOT_2 = 0;
	double max_Vx = 0;
	double max_Vy = 0;
	double max_Vz = 0;
	int max_weight_Vx = 0;
	int max_weight_Vy = 0;
	int max_weight_Vz = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz / static_cast<double>(std::pow(2.0, weights[2][i])));
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_e_V_TOT_2 = e_V_TOT_2;
			max_V_TOT_2 = V_TOT_2;
			max_Vx = Vx[i];
			max_Vy = Vy[i];
			max_Vz = Vz[i];
			max_index = i;
			max_weight_Vx = weights[0][i];
			max_weight_Vy = weights[1][i];
			max_weight_Vz = weights[2][i];
			// max_weight_index = weight_index;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		// estimate
		auto i = max_index;
		double estimate_error = max_value;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5;		        		
			double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz / static_cast<double>(std::pow(2.0, weights[2][i])));
			// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			// if(ebs[0]/eb_Vx > 5) break;
			if(ebs[0] / eb_Vx > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	
	return true;
}

template<class T>
bool halfing_error_V_TOT_coordinate(const T * Vx, const T * Vy, const T * Vz, size_t n, const double tau, std::vector<double>& ebs, std::vector<std::vector<int>> weights){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	double max_e_V_TOT_2 = 0;
	double max_V_TOT_2 = 0;
	double max_Vx = 0;
	double max_Vy = 0;
	double max_Vz = 0;
	int max_weight_Vx = 0;
	int max_weight_Vy = 0;
	int max_weight_Vz = 0;
	// int weight_index = 0;
	// int max_weight_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz / static_cast<double>(std::pow(2.0, weights[2][i])));
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_e_V_TOT_2 = e_V_TOT_2;
			max_V_TOT_2 = V_TOT_2;
			max_Vx = Vx[i];
			max_Vy = Vy[i];
			max_Vz = Vz[i];
			max_index = i;
			max_weight_Vx = weights[0][i];
			max_weight_Vy = weights[1][i];
			max_weight_Vz = weights[2][i];
			// max_weight_index = weight_index;
		}
	}
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_V_TOT_2 = " << max_e_V_TOT_2 << ", VTOT_2 = " << max_V_TOT_2 << ", Vx = " << max_Vx << ", Vy = " << max_Vy << ", Vz = " << max_Vz << std::endl;
	// std::cout << "max_weight_Vx = " << max_weight_Vx << ", max_weight_Vy = " << max_weight_Vy << ", max_weight_Vz = " << max_weight_Vz << std::endl;

	// ************************************************ Vx, Vy and Vz ************************************************
	if(max_value > tau){
		// estimate
		auto i = max_index;
		double estimate_error = max_value;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		while(estimate_error > tau){
			// std::cout << "coordinate decrease\n";
			// std::cout << "coordinate decrease, eb_Vx / ebs[0] = " << eb_Vx / ebs[0] << std::endl;
    		double estimate_error_Vx = 0;
			{
				double eb_Vx_ = eb_Vx / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx_ / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz / static_cast<double>(std::pow(2.0, weights[2][i])));
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vx = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			double estimate_error_Vy = 0;
			{
				double eb_Vy_ = eb_Vy / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy_ / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz / static_cast<double>(std::pow(2.0, weights[2][i])));
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vy = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			double estimate_error_Vz = 0;
			{
				double eb_Vz_ = eb_Vz / 1.5;
				double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx / static_cast<double>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square((double) Vy[i], eb_Vy / static_cast<double>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square((double) Vz[i], eb_Vz_ / static_cast<double>(std::pow(2.0, weights[2][i])));
				// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
				estimate_error_Vz = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			}
			// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << std::endl;
			const double relative_epsilon = 1e-3;
			double min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz});
			double epsilon = std::max(relative_epsilon * min_error, static_cast<double>(1e-12));
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
			if (ebs[0] / eb_Vx > 10 || ebs[1] / eb_Vy > 10 || ebs[2] / eb_Vz > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
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
std::vector<size_t> retrieve_V_TOT_Dummy(std::string rdata_file_prefix, double tau, std::vector<double> ebs, size_t num_elements, int weighted, double & max_act_error, double & max_est_error, size_t & weight_file_size, bool decrease_method){
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
			std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
			std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			MDR::print_vec(ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
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
			reconstructors.back().load_metadata();
			reconstructors.back().load_weight();
			reconstructors.back().span_weight();
			weights[i] = reconstructors.back().get_int_weights();
		}
		weight_file_size = reconstructors[0].get_weight_file_size();
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
			std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
			std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			MDR::print_vec(ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
	}
	return total_retrieved_size;
}

std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> sz3_reconstructors;

template<class T>
std::vector<size_t> retrieve_V_TOT_SZ3(std::string rdata_file_prefix, double tau, std::vector<double> ebs, size_t num_elements, int weighted, double & max_act_error, double & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
	if(!weighted){
		std::vector<PDR::ApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::NegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
		for(int i=0; i<n_variable; i++){
			std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
			else tolerance_met = halfing_error_V_TOT_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
        local_elapsed_time += MPI_Wtime();
	}
	else{
		std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::SZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
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
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
        local_elapsed_time += MPI_Wtime();
		sz3_reconstructors = reconstructors;
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT_PMGARD(std::string rdata_file_prefix, double tau, std::vector<double> ebs, size_t num_elements, int weighted, double & max_act_error, double & max_est_error, size_t & weight_file_size, bool decrease_method){
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
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
			std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			MDR::print_vec(ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
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
		weight_file_size = reconstructors[0].get_weight_file_size();
		while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), -1);
				// auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
				total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			for(int i=0; i<n_variable; i++){
				std::cout << total_retrieved_size[i] << ", ";
			}
			std::cout << ": total_size = " << std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) << std::endl;

			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
			MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
			MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			MDR::print_vec(ebs);
			tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
			std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
			MDR::print_vec(ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
	}
	return total_retrieved_size;
}

template<class T>
std::vector<size_t> retrieve_V_TOT_GE(std::string rdata_file_prefix, double tau, std::vector<double> ebs, size_t num_elements, int weighted, double & max_act_error, double & max_est_error, size_t & weight_file_size, bool decrease_method){
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
			reconstructors.back().load_metadata();
    	}
        local_elapsed_time = -MPI_Wtime();
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
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
			else tolerance_met = halfing_error_V_TOT_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
        local_elapsed_time += MPI_Wtime();
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
				// auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
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
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
			// MDR::print_vec(ebs);
			if(!decrease_method) tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
			else tolerance_met = halfing_error_V_TOT_coordinate(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs, weights);
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
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
			max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
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

	using T = float;
	using T_stream = uint32_t;
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
	int exp = static_cast<int>(std::round(std::log10(target_rel_eb)));
	std::string wdata_file_prefix = output_path + "/1e" + std::to_string(exp) + "/";

    size_t num_elements = 0;
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);

    std::vector<double> ebs;
	ebs.push_back(compute_global_value_range(Vx_ori)*target_rel_eb);
	ebs.push_back(compute_global_value_range(Vy_ori)*target_rel_eb);
	ebs.push_back(compute_global_value_range(Vz_ori)*target_rel_eb);
	// ebs.push_back(compute_max_abs_value(Vx_ori.data(), Vx_ori.size()*target_rel_eb));
	// ebs.push_back(compute_max_abs_value(Vy_ori.data(), Vy_ori.size()*target_rel_eb));
	// ebs.push_back(compute_max_abs_value(Vz_ori.data(), Vz_ori.size()*target_rel_eb));
	int n_variable = ebs.size();

    std::vector<T> V_TOT(num_elements);
    compute_VTOT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT.data());

    double tau, global_tau;
	global_tau = compute_global_value_range(V_TOT) * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
	V_TOT_ori = V_TOT.data();

	double local_max_act_error = 0;
	double local_max_est_error = 0;
	size_t weight_file_size = 0;
	std::vector<size_t> total_retrieved_size;
	switch (compressor)
	{
	case Dummy:
		total_retrieved_size = retrieve_V_TOT_Dummy<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	case SZ3:
		total_retrieved_size = retrieve_V_TOT_SZ3<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	case PMGARD:
		total_retrieved_size = retrieve_V_TOT_PMGARD<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
		break;
	case GE:
		total_retrieved_size = retrieve_V_TOT_GE<T>(rdata_file_prefix, tau, ebs, num_elements, weighted, local_max_act_error, local_max_est_error, weight_file_size, decrease_method);
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
	// 	if(rank == i) {
	// 		std::cout << "rank = " << rank << " act_iter = " << iter << ", local_bitrate = " << 8 * local_total_size * 1.0 / (num_elements * n_variable) << std::endl;
	// 	}
	// }

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);

	size_t total_2 = 0;
    for(int i=0; i<n_variable; i++){
        auto tmp_count = sz3_reconstructors[i].get_offsets();
		auto approximator_count = sz3_reconstructors[i].get_approximator_size();
		auto metadata_count = sz3_reconstructors[i].get_metadata_size();

		std::vector<unsigned long long int> count(tmp_count.begin(), tmp_count.end());
        auto offsets(count);
        for(int j=0; j<offsets.size(); j++) offsets[j] = 0;
		unsigned long long int approximator_offset = 0;
		unsigned long long int metadata_offset = 0;

        auto buffer(offsets);
		unsigned long long int approximator_buffer;
		unsigned long long int metadata_buffer;

        for(int j=0; j<size; j++){
            if(j == rank){
                if(j != 0) {
                    MPI_Recv(&offsets[0], offsets.size(), MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&approximator_offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&metadata_offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                for(int k=0; k<offsets.size(); k++){
                    buffer[k] = offsets[k] + count[k];
                }
				approximator_buffer = approximator_offset + approximator_count;
				metadata_buffer = metadata_offset + metadata_count;

                if(j != size - 1) {
					MPI_Send(&buffer[0], offsets.size(), MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
					MPI_Send(&approximator_buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
					MPI_Send(&metadata_buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
				}
            }
        }
		std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";
        for(int k=0; k<offsets.size(); k++){
            std::string file_level = rdir_prefix + "level_" + std::to_string(k) + ".bin";
            size_t num_char = 0;
            auto level_data = MGARD::readfile<unsigned char>(file_level.c_str(), num_char);
            MPI_File file;
            std::string filename = wdata_file_prefix + var_name_out[i] + "_aggregated_level_" + std::to_string(k) + ".dat";
            MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
            MPI_File_write_at(file, offsets[k], level_data.data(), count[k], MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
            MPI_File_close(&file);
            total_2 += count[k];
        }
		std::string approximator_path = rdir_prefix + "approximator.dat";
		std::string metadata_path = rdir_prefix + "metadata.bin";
		size_t approximator_num_char = 0;
		size_t metadata_num_char = 0;
		auto approximator_data = MGARD::readfile<char>(approximator_path.c_str(), approximator_num_char);
		auto metadata_data = MGARD::readfile<unsigned char>(metadata_path.c_str(), metadata_num_char);
		MPI_File approximator_file;
		MPI_File metadata_file;
		std::string approximator_filename = wdata_file_prefix + var_name_out[i] + "_aggregated_approximator.dat";
		std::string metadata_filename = wdata_file_prefix + var_name_out[i] + "_aggregated_metadata.bin";
		MPI_File_open(MPI_COMM_WORLD, approximator_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &approximator_file);
		MPI_File_write_at(approximator_file, approximator_offset, approximator_data.data(), approximator_count, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&approximator_file);
		MPI_File_open(MPI_COMM_WORLD, metadata_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &metadata_file);
		MPI_File_write_at(metadata_file, metadata_offset, metadata_data.data(), metadata_count, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&metadata_file);
		if(i == 0){
			auto weight_count = sz3_reconstructors[i].get_weight_file_size();
			unsigned long long int weight_offset = 0;
			unsigned long long int weight_buffer; 
			for(int j=0; j<size; j++){
				if(j == rank){
					if(j != 0) {
						MPI_Recv(&weight_offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
					weight_buffer = weight_offset + weight_count;
					if(j != size - 1) {
						MPI_Send(&weight_buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
					}
				}
			}
			std::string weight_path = rdir_prefix + "weight.bin";
			size_t weight_num_char = 0;
			auto weight_data = MGARD::readfile<unsigned char> (weight_path.c_str(), weight_num_char);
			MPI_File weight_file;
			std::string weight_filename = wdata_file_prefix + var_name_out[i] + "_aggregated_weight.bin";
			MPI_File_open(MPI_COMM_WORLD, weight_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &weight_file);
			MPI_File_write_at(weight_file, weight_offset, weight_data.data(), weight_count, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
			MPI_File_close(&weight_file);
		}
    }

    MPI_Finalize();

    return 0;
}