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

using namespace MDR;

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
float * Temp_ori = NULL;
std::vector<float> error_Temp;
std::vector<float> error_est_Temp;


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
	std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << ", e_T = " << max_e_T << ", T = " << max_T << ", P = " << max_P << ", D = " << max_D << std::endl;
    std::cout << "max_weight_P = " << max_weight_P << ", max_weight_D = " << max_weight_D << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		T estimate_error = max_value;
        // float T = c_1 * P[i] / D[i];
        T eb_P = ebs[0];
        T eb_D = ebs[1];
        while(estimate_error > tau){
    		std::cout << "uniform decrease\n";
            eb_P = eb_P / 1.5;
            eb_D = eb_D / 1.5;
            estimate_error = c_1 * compute_bound_division(P[i], D[i], eb_P / static_cast<T>(std::pow(2.0, weights[0][i])), eb_D / static_cast<T>(std::pow(2.0, weights[1][i])));
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}
	return true;
}

int main(int argc, char ** argv){

    using T = float;
	using T_stream = typename std::conditional<std::is_same<T, double>::value, uint64_t, uint32_t>::type;
	int argv_id = 1;
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

    std::vector<T> Temp(num_elements);
    compute_T(P_ori.data(), D_ori.data(), num_elements, Temp.data());
	Temp_ori = Temp.data();
    T tau = compute_value_range(Temp)*target_rel_eb;
    // /* test
    std::string filename = "./Result/Temp_ori.dat";
    std::ofstream outfile1(filename, std::ios::binary);
    if (!outfile1.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
		exit(-1);
    }
	std::cout << "Temp.size: " << Temp.size() << std::endl;
    outfile1.write(reinterpret_cast<const char*>(Temp.data()), Temp.size() * sizeof(float));
    
    outfile1.close();
    std::cout << "Data saved successfully to " << filename << std::endl;
    //*/

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
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<float>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);

    int iter = 0;
    int max_iter = 10;
	bool tolerance_met = false;
	double max_act_error = 0, max_est_error = 0;
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
	    error_Temp = std::vector<float>(num_elements);
	    error_est_Temp = std::vector<float>(num_elements);
		std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    MDR::print_vec(ebs);
	    tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs, weights);
		std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    MDR::print_vec(ebs);
        // /* test
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
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
	std::cout << "iter = " << iter << std::endl;
   
   	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0);
    total_size = total_size + reconstructors[0].get_weight_file_size();
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_size[i] << ", ";
    }
    std::cout << "weight_file_size: " << reconstructors[0].get_weight_file_size();
    std::cout << std::endl;
	// MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}