#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <mpi.h>
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

const std::string data_prefix_path = "/pscratch/xli281_uksr/xliang/GE";

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
double * PT_ori = NULL;
std::vector<double> error_PT;
std::vector<double> error_est_PT;
int iter = 0;

template <class T>
T print_max_abs(int rank, const std::string& name, const std::vector<T>& vec){
	T max = fabs(vec[0]);
	for(int i=1; i<vec.size(); i++){
		if(max < fabs(vec[i])) max = fabs(vec[i]);
	}
	// printf("Processor %d var %s: max absolute value =  %.4f\n", rank, name.c_str(), max);
	return max;
}

template<class T>
bool halfing_error_PT_uniform(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const T tau, std::vector<T>& ebs, std::vector<std::vector<int>> weights){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	T eb_P = ebs[3];
	T eb_D = ebs[4];
	T R = 287.1;
	T gamma = 1.4;
	T mi = 3.5;
	T mu_r = 1.716e-5;
	T T_r = 273.15;
	T S = 110.4;
	T c_1 = 1.0 / R;
	T c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	T max_value = 0;
	int max_index = 0;
    T max_Vx = 0;
    T max_Vy = 0;
    T max_Vz = 0;
    T max_P = 0;
    T max_D = 0;
    T max_e_VTOT_2 = 0;
    T max_VTOT_2 = 0;
    T max_e_T = 0;
    T max_T = 0;
    T max_C = 0;
    T max_e_C = 0;
    T max_e_Mach = 0;
    T max_Mach = 0;
    T max_e_Mach_tmp_mi = 0;
    T max_Mach_tmp_mi = 0;
    T max_e_PT = 0;
    T max_PT = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
		T Temp = P[i] / (D[i] * R);
		T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		T C = c_2 * sqrt(Temp);
		T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		T Mach = V_TOT / C;
		T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
		T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
		T e_Mach_tmp_mi = 0;
		for(int i=1; i<=7; i++){
			e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
		}
		T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
		T e_PT = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
		T PT = P[i] * Mach_tmp_mi;

		error_est_PT[i] = e_PT;
		error_PT[i] = PT - PT_ori[i];
		if(max_value < error_est_PT[i]){
			max_value = error_est_PT[i];
			max_index = i;
            max_Vx = Vx[i];
            max_Vy = Vy[i];
            max_Vz = Vz[i];
            max_P = P[i];
            max_D = D[i];
            max_e_VTOT_2 = e_V_TOT_2;
            max_VTOT_2 = V_TOT_2;
            max_e_T = e_T;
            max_T = Temp;
            max_e_C = e_C;
            max_C = C;
            max_e_Mach = e_Mach;
            max_Mach = Mach;
            max_e_Mach_tmp_mi = e_Mach_tmp_mi;
            max_Mach_tmp_mi = Mach_tmp_mi;
            max_e_PT = e_PT;
            max_PT = PT;
		}
	}
    if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
		T eb_Vx = ebs[0];
		T eb_Vy = ebs[1];
		T eb_Vz = ebs[2];
		T eb_P = ebs[3];
		T eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            T estimate_error_Vx = 0;
            {
                T eb_Vx_ = eb_Vx / 1.5;
                T e_V_TOT_2 = 0;
                if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx_ / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
                T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
                T e_V_TOT = 0;
                if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
                T V_TOT = sqrt(V_TOT_2);
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
                T Temp = P[i] / (D[i] * R);
                T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
                T C = c_2 * sqrt(Temp);
                T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
                T Mach = V_TOT / C;
                T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
                T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
                T e_Mach_tmp_mi = 0;
                for(int i=1; i<=7; i++){
                    e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
                }
                T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
                estimate_error_Vx = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
            }
            T estimate_error_Vy = 0;
            {
                T eb_Vy_ = eb_Vy / 1.5;
                T e_V_TOT_2 = 0;
                if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy_ / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
                T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
                T e_V_TOT = 0;
                if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
                T V_TOT = sqrt(V_TOT_2);
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
                T Temp = P[i] / (D[i] * R);
                T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
                T C = c_2 * sqrt(Temp);
                T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
                T Mach = V_TOT / C;
                T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
                T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
                T e_Mach_tmp_mi = 0;
                for(int i=1; i<=7; i++){
                    e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
                }
                T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
                estimate_error_Vy = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
            }
            T estimate_error_Vz = 0;
            {
                T eb_Vz_ = eb_Vz / 1.5;
                T e_V_TOT_2 = 0;
                if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz_ / static_cast<T>(std::pow(2.0, weights[2][i])));
                T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
                T e_V_TOT = 0;
                if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
                T V_TOT = sqrt(V_TOT_2);
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
                T Temp = P[i] / (D[i] * R);
                T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
                T C = c_2 * sqrt(Temp);
                T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
                T Mach = V_TOT / C;
                T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
                T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
                T e_Mach_tmp_mi = 0;
                for(int i=1; i<=7; i++){
                    e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
                }
                T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
                estimate_error_Vz = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
            }
            T estimate_error_P = 0;
            {
                T eb_P_ = eb_P / 1.5;
                T e_V_TOT_2 = 0;
                if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
                T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
                T e_V_TOT = 0;
                if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
                T V_TOT = sqrt(V_TOT_2);
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P_ / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D / static_cast<T>(std::pow(2.0, weights[4][i])));
                T Temp = P[i] / (D[i] * R);
                T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
                T C = c_2 * sqrt(Temp);
                T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
                T Mach = V_TOT / C;
                T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
                T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
                T e_Mach_tmp_mi = 0;
                for(int i=1; i<=7; i++){
                    e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
                }
                T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
                estimate_error_P = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P_ / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
            }
            T estimate_error_D = 0;
            {
                T eb_D_ = eb_D / 1.5;
                T e_V_TOT_2 = 0;
                if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx / static_cast<T>(std::pow(2.0, weights[0][i]))) + compute_bound_x_square(Vy[i], eb_Vy / static_cast<T>(std::pow(2.0, weights[1][i]))) + compute_bound_x_square(Vz[i], eb_Vz / static_cast<T>(std::pow(2.0, weights[2][i])));
                T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
                T e_V_TOT = 0;
                if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
                T V_TOT = sqrt(V_TOT_2);
                T e_T = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), eb_D_ / static_cast<T>(std::pow(2.0, weights[4][i])));
                T Temp = P[i] / (D[i] * R);
                T e_C = c_2*compute_bound_square_root_x(Temp, e_T);
                T C = c_2 * sqrt(Temp);
                T e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
                T Mach = V_TOT / C;
                T e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
                T Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
                T e_Mach_tmp_mi = 0;
                for(int i=1; i<=7; i++){
                    e_Mach_tmp_mi += C7i[i] * pow(Mach_tmp, 7-i) * pow(e_Mach_tmp, i);
                }
                T Mach_tmp_mi = sqrt(pow(Mach_tmp, 7));
                estimate_error_D = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P / static_cast<T>(std::pow(2.0, weights[3][i])), e_Mach_tmp_mi);
            }
    		// std::cout << estimate_error_Vx << " " << estimate_error_Vy << " " << estimate_error_Vz << " " << estimate_error_P << " " << estimate_error_D << std::endl;
			const T relative_epsilon = 1e-3;
            T min_error = std::min({estimate_error_Vx, estimate_error_Vy, estimate_error_Vz, estimate_error_P, estimate_error_D});
            T epsilon = std::max(relative_epsilon * min_error, static_cast<T>(1e-12));
            bool close_Vx = fabs(estimate_error_Vx - min_error) < epsilon;
            bool close_Vy = fabs(estimate_error_Vy - min_error) < epsilon;
            bool close_Vz = fabs(estimate_error_Vz - min_error) < epsilon;
            bool close_P  = fabs(estimate_error_P - min_error) < epsilon;
            bool close_D  = fabs(estimate_error_D - min_error) < epsilon;
            estimate_error = min_error;
            if (close_Vx) eb_Vx /= 1.5;
            if (close_Vy) eb_Vy /= 1.5;
            if (close_Vz) eb_Vz /= 1.5;
            if (close_P)  eb_P /= 1.5;
            if (close_D)  eb_D /= 1.5;
            if ((ebs[0] / eb_Vx > 10) || (ebs[1] / eb_Vy > 10) || (ebs[2] / eb_Vz > 10) || (ebs[3] / eb_P > 10) || (ebs[4] / eb_D > 10)) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
		return false;
	}

	return true;
}

template<class T>
std::vector<size_t> retrieve_PT_GE(int rank, T tau, std::vector<T> ebs, size_t num_elements, const std::vector<unsigned char>& mask, T & max_act_error, T & max_est_error, size_t & weight_file_size){
    int max_iter = 30;
    bool tolerance_met = false;
    int n_variable = ebs.size();
    // std::cout << "n_variable = " << n_variable << std::endl; 
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);
    std::string block_path = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/block_sizes.dat";

    std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + varlist[i];
        std::string metadata_file = rdir_prefix + "/metadata.bin";
        std::string weight_file = rdir_prefix + "/weight.bin";
        std::vector<std::string> files;
        int num_levels = 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::GEApproximator<T>();
        auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
        auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<T>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        reconstructors.push_back(generateWBPReconstructor_GE<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
        if(i < 3) reconstructors.back().mask = mask;
        reconstructors.back().load_metadata();
        reconstructors.back().load_weight(weight_file);
        reconstructors.back().span_weight(block_path);
        weights[i] = reconstructors.back().get_int_weights();
    }
    weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[n_variable - 1].get_weight_file_size();
    // std::cout << "rank = " << rank << ", weight_file_size = " << weight_file_size << std::endl;
    while((!tolerance_met) && (iter < max_iter)){
        iter ++;
        for(int i=0; i<n_variable; i++){
            std::string approximator_path = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + varlist[i] + "/approximator.dat";
            auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), block_path, approximator_path, -1);
            total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
        }
        Vx_dec = reconstructed_vars[0].data();
        Vy_dec = reconstructed_vars[1].data();
        Vz_dec = reconstructed_vars[2].data();
        P_dec = reconstructed_vars[3].data();
        D_dec = reconstructed_vars[4].data();
        // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
        // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
        // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
        // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
        // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
        error_PT = std::vector<T>(num_elements);
        error_est_PT = std::vector<T>(num_elements);
        // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
        // MDR::print_vec(ebs);
        tolerance_met = halfing_error_PT_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
        // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
        // MDR::print_vec(ebs);
        max_act_error = print_max_abs(names[1] + " error", error_PT);
        max_est_error = print_max_abs(names[1] + " error_est", error_est_PT);   	
    }

    return total_retrieved_size;
}

int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    using T = double;
    using T_stream = uint64_t;
    T target_rel_eb = atof(argv[1]);
    const int target_level = 8;

    size_t num_elements = 0;
	std::string filename = "/pscratch/xli281_uksr/xliang/GE/sol_4114800_aver_b" + std::to_string(rank) + ".bp/"; 
	Vx_ori = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityX.dat.rod").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityY.dat.rod").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/VelocityZ.dat.rod").c_str(), num_elements);
    P_ori = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/Pressure.dat.rod").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_prefix_path + "/sol_4114800_aver_b" + std::to_string(rank) + ".bp/Density.dat.rod").c_str(), num_elements);

    std::vector<T> ebs;
    ebs.push_back(compute_value_range(Vx_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vy_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vz_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(D_ori)*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Vx_ori, Vy_ori, Vz_ori, P_ori, D_ori};
    std::vector<T> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

	struct timespec start, end;
	int err;
	T elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> PT(num_elements);
    compute_PT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), P_ori.data(), D_ori.data(), num_elements, PT.data());
	PT_ori = PT.data();
    // T tau = compute_value_range(PT)*target_rel_eb;

    double tau;
	double local_max = -9999, local_min = 9999;
	double global_max = 0, global_min = 0;
	for(int i=0; i<num_elements; i++){
		if(PT[i] > local_max) local_max = PT[i];
		if(PT[i] < local_min) local_min = PT[i];
	}
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	double max_value_range = global_max - global_min;
	tau = max_value_range*target_rel_eb;

    std::string block_path = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/block_sizes.dat";
    std::string mask_file = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/mask.bin";
    uint32_t mask_file_size = 0;
    auto mask = readmask(mask_file.c_str(), mask_file_size);
    T max_act_error = 0, max_est_error = 0;
    size_t weight_file_size = 0;
    // std::vector<size_t> total_retrieved_size;

    // total_retrieved_size = retrieve_PT_GE(rank, tau, ebs, num_elements, mask, max_act_error, max_est_error, weight_file_size);

    int max_iter = 30;
    bool tolerance_met = false;
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
    std::vector<size_t> total_retrieved_size(n_variable, 0);

    std::vector<PDR::GEReconstructor<T, PDR::GEApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, T_stream>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + varlist[i];
        std::string metadata_file = rdir_prefix + "/metadata.bin";
        std::string weight_file = rdir_prefix + "/weight.bin";
        std::vector<std::string> files;
        int num_levels = 1;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto approximator = PDR::GEApproximator<T>();
        auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
        auto compressor = AdaptiveLevelCompressor(32);
        // auto encoder = WeightedNegaBinaryBPEncoder<T, T_stream>();
        // auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<T>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        reconstructors.push_back(generateWBPReconstructor_GE<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
        if(i < 3) reconstructors.back().mask = mask;
        reconstructors.back().load_metadata();
        reconstructors.back().load_weight(weight_file);
        reconstructors.back().span_weight(block_path);
        weights[i] = reconstructors.back().get_int_weights();
    }
    weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[n_variable - 1].get_weight_file_size();
    while((!tolerance_met) && (iter < max_iter)){
        iter ++;
        for(int i=0; i<n_variable; i++){
            std::string approximator_path = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + varlist[i] + "/approximator.dat";
            auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i] / static_cast<T>(std::pow(2.0, reconstructors[i].get_max_weight())), block_path, approximator_path, -1);
            total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
        }
        Vx_dec = reconstructed_vars[0].data();
        Vy_dec = reconstructed_vars[1].data();
        Vz_dec = reconstructed_vars[2].data();
        P_dec = reconstructed_vars[3].data();
        D_dec = reconstructed_vars[4].data();
        error_PT = std::vector<T>(num_elements);
        error_est_PT = std::vector<T>(num_elements);
        tolerance_met = halfing_error_PT_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs, weights);
        max_act_error = print_max_abs(names[1] + " error", error_PT);
        max_est_error = print_max_abs(names[1] + " error_est", error_est_PT);   	
    }

	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (T)(end.tv_sec - start.tv_sec) + (T)(end.tv_nsec - start.tv_nsec)/(T)1000000000;

   	size_t total_size = mask_file_size + std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0) + weight_file_size;
	T cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;

	int maxit;
	MPI_Reduce(&iter, &maxit, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) std::cout << "overall max iter = " << maxit << std::endl;
	// std::cout << "rank " << rank << " max iter = " << iter << std::endl;
    double max_time;
	elapsed_time += MPI_Wtime();
	MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // printf("Processor %d total size = %d\n", rank, total_size);
    // MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) printf("Target PT error = %.4f\n", tau);
	double max_error = 0;
	max_error = print_max_abs(rank, "PT error", error_PT);
	double max_pt_error = 0;
	MPI_Reduce(&max_error, &max_pt_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated PT error = %.4f\n", max_pt_error);
	max_error = print_max_abs(rank, "PT error", error_est_PT);
	double max_pt_error_est = 0;
	MPI_Reduce(&max_error, &max_pt_error_est, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated PT est error = %.4f\n", max_pt_error_est);
	unsigned long long int total_num = 0;
	MPI_Reduce(&num_elements, &total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int total_retrieved = 0;
	MPI_Reduce(&total_size, &total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.4f, retrieved_size = %ld, total_num_elements = %ld\n", 8*total_retrieved * 1.0 / (total_num * n_variable), total_retrieved, total_num);
	if(!rank) printf("elapsed_time = %.6f\n", max_time);

	size_t total_2 = 0;
    for(int i=0; i<n_variable; i++){
    	auto count = reconstructors[i].get_offsets();
		auto offsets(count);
		auto buffer(count);
		for(int j=0; j<offsets.size(); j++) offsets[j] = 0;
		for(int j=0; j<size; j++){
			if(j == rank){
				if(j != 0) {
					MPI_Recv(&buffer[0], offsets.size(), MPI_UNSIGNED, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(int k=0; k<offsets.size(); k++){
						offsets[k] = buffer[k] + count[k];
					}
				}
				if(j != size - 1) MPI_Send(&offsets[0], offsets.size(), MPI_UNSIGNED, j+1, 0, MPI_COMM_WORLD);
			}
		}
		for(int k=0; k<offsets.size(); k++){
			std::string rdir_prefix = "/pscratch/xli281_uksr/xliang/GE/block_" + std::to_string(rank) + "_refactored/" + varlist[i] + "/";
			// printf("Processor %d offset = %d, count = %d, offset + count = %d\n", rank, offsets[k], count[k], offsets[k] + count[k]);
			std::string file_level = rdir_prefix + "level_" + std::to_string(k) + ".bin";
			size_t num_char = 0;
			auto level_data = MGARD::readfile<unsigned char>(file_level.c_str(), num_char);
			MPI_File file;
            std::string filename = "/pscratch/xli281_uksr/xliang/GE/PT_QPRO_" + std::to_string(target_rel_eb) + "_" + varlist[i] + "_aggregated_level_" + std::to_string(k) + ".dat";
			MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
			MPI_File_write_at(file, offsets[k], level_data.data(), count[k], MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
			MPI_File_close(&file);
			total_2 += count[k];
		}
    }

	// std::cout << "requested error = " << tau << std::endl;
	// std::cout << "max_est_error = " << max_est_error << std::endl;
	// std::cout << "max_act_error = " << max_act_error << std::endl;
	// std::cout << "iter = " << iter << std::endl;
   
	// std::cout << "each retrieved size:";
    // for(int i=0; i<n_variable; i++){
    //     std::cout << total_retrieved_size[i] << ", ";
    // }
	// MDR::print_vec(total_retrieved_size);
	// std::cout << "aggregated cr = " << cr << std::endl;
	// printf("elapsed_time = %.6f\n", elapsed_time);

    MPI_Finalize();

    return 0;
}