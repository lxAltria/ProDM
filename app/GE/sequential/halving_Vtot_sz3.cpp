#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "MDR/Synthesizer4GE.hpp"

using namespace MDR;
using namespace QoI;

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
double * V_TOT_ori = NULL;
std::vector<double> error_V_TOT;
std::vector<double> error_est_V_TOT;


template<class T>
bool halfing_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const std::vector<unsigned char>& mask, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_index = i;
		}

	}
	std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
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
    		std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 							        		
			double e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
			estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	return true;
}

int main(int argc, char** argv){
    using T = double;
	int argv_id = 1;
    double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Vx_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vy_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vz_ori)*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Vx_ori, Vy_ori, Vz_ori};
    std::vector<double> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> V_TOT(num_elements);
    compute_VTOT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT.data());
	V_TOT_ori = V_TOT.data();
    double tau = compute_value_range(V_TOT)*target_rel_eb;

    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);

    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	std::vector<size_t> total_retrieved_sizes(n_variable, 0);

    int iter = 0;
    int max_iter = 5;
	bool tolerance_met = false;
	double max_act_error = 0, max_est_error = 0;
    std::vector<int> current_ind(n_variable, -1);
    T * reconstructed_data = (T *) malloc(num_elements * sizeof(T));
    while((!tolerance_met) && (iter < max_iter)){
        iter ++;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i];
            double file_eb = 0.1;
            auto file_ind = find_index(ebs[i]/var_range[i], file_eb);
            std::cout << "file_ind = " << file_ind << std::endl;
            std::cout << "Requested relative tolerance = " << ebs[i]/var_range[i] << ", expected relative tolerance = " << file_eb << "\n"; 
			// if(file_ind >= 8){
			// 	std::cout << "reached maximum capacity" << std::endl;
			// 	file_ind = 7;
			// }			
            std::string filename = rdir_prefix + "_refactored/SZ3_eb_" + std::to_string(file_ind) + ".bin";
            size_t n = 0;
            auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
            SZ3_decompress(cmpData.data(), n, reconstructed_data);
            if(file_ind > current_ind[i]){
                current_ind[i] = file_ind;
                total_retrieved_sizes[i] += n;
            }
			// reconstruct with mask
			int index = 0;
			for(int j=0; j<num_elements; j++){
				if(mask[j]){
					reconstructed_vars[i][j] = reconstructed_data[index ++];
				}
				else reconstructed_vars[i][j] = 0;
			}
			std::cout << varlist[i] << " bitrate = " <<  total_retrieved_sizes[i] * 1.0 / num_valid_data << std::endl;
        }        
        Vx_dec = reconstructed_vars[0].data();
        Vy_dec = reconstructed_vars[1].data();
        Vz_dec = reconstructed_vars[2].data();
	    // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
	    // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
	    // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
	    error_V_TOT = std::vector<double>(num_elements);
	    error_est_V_TOT = std::vector<double>(num_elements);
		std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, mask, tau, ebs);
		std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    // std::cout << names[0] << " requested error = " << tau << std::endl;
	    max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT); 
	    max_act_error = print_max_abs(names[0] + " actual error", error_V_TOT);
    }
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
    free(reconstructed_data);
	std::cout << "iter = " << iter << std::endl;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_sizes[i] << ", ";
    }
    std::cout << std::endl;
	// MDR::print_vec(total_retrieved_size);
	size_t total_size = std::accumulate(total_retrieved_sizes.begin(), total_retrieved_sizes.end(), 0);
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "aggregated cr = " << cr << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}