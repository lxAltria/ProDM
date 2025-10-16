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
#define HPEZ 4
using namespace MDR;

const std::vector<std::string> species = {"H2", "O2", "H2O", "H", "O", "OH", "HO2", "H2O2"};

std::vector<double> Xi_ori;
std::vector<double> Xj_ori;
double * Xi_dec = NULL;
double * Xj_dec = NULL;
double * XiXj_ori = NULL;
std::vector<double> error_XiXj;
std::vector<double> error_est_XiXj;
int iter = 0;

template<class T>
bool halfing_error_XiXj_uniform(const T * Xi, const T * Xj, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Xi = ebs[0];
	double eb_Xj = ebs[1];
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_XiXj = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
        double XiXj = Xi[i] * Xj[i];

		error_est_XiXj[i] = e_XiXj;
		error_XiXj[i] = XiXj - XiXj_ori[i];
		if(max_value < error_est_XiXj[i]){
			max_value = error_est_XiXj[i];
			max_index = i;
		}
	}
	// std::cout << "XiXj: max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Xi = ebs[0];
		double eb_Xj = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Xi = eb_Xi / 1.5;
			eb_Xj = eb_Xj / 1.5;
			estimate_error = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
		}
		ebs[0] = eb_Xi;
		ebs[1] = eb_Xj;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_XiXj_coordinate(const T * Xi, const T * Xj, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Xi = ebs[0];
	double eb_Xj = ebs[1];
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_XiXj = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
        double XiXj = Xi[i] * Xj[i];

		error_est_XiXj[i] = e_XiXj;
		error_XiXj[i] = XiXj - XiXj_ori[i];
		if(max_value < error_est_XiXj[i]){
			max_value = error_est_XiXj[i];
			max_index = i;
		}
	}
	// std::cout << "XiXj: max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Xi = ebs[0];
		double eb_Xj = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            double estimate_error_Xi = 0;
            {
                double eb_Xi_ = eb_Xi / 1.5;
                estimate_error_Xi = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi_, eb_Xj);
            }
            double estimate_error_Xj = 0;
            {
                double eb_Xj_ = eb_Xj / 1.5;
                estimate_error_Xj = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj_);
            }
			double sum_err = 2 * estimate_error - (estimate_error_Xi + estimate_error_Xj);
            double w_Xi = (estimate_error - estimate_error_Xi) / sum_err;
            double w_Xj = (estimate_error - estimate_error_Xj) / sum_err;
            // === Smooth proportional update ===
            double factor_base = 1.5;
            double alpha = 1.0 - (1.0 / factor_base);
            eb_Xi *= (1.0 - alpha * w_Xi);
            eb_Xj *= (1.0 - alpha * w_Xj);
			estimate_error = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
		}
		ebs[0] = eb_Xi;
		ebs[1] = eb_Xj;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_XiXj_uniform(const T * Xi, const T * Xj, size_t n, const double tau, std::vector<double>& ebs, std::vector<std::vector<int>>& weights){
	double eb_Xi = ebs[0];
	double eb_Xj = ebs[1];
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_XiXj = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi, -weights[0][i]), ldexp(eb_Xj, -weights[1][i]));
        double XiXj = Xi[i] * Xj[i];

		error_est_XiXj[i] = e_XiXj;
		error_XiXj[i] = XiXj - XiXj_ori[i];
		if(max_value < error_est_XiXj[i]){
			max_value = error_est_XiXj[i];
			max_index = i;
		}
	}
	// std::cout << "XiXj: max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Xi = ebs[0];
		double eb_Xj = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Xi = eb_Xi / 1.5;
			eb_Xj = eb_Xj / 1.5;
			estimate_error = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi, -weights[0][i]), ldexp(eb_Xj, -weights[1][i]));
		}
		ebs[0] = eb_Xi;
		ebs[1] = eb_Xj;
		return false;
	}
	return true;
}

template<class T>
bool halfing_error_XiXj_coordinate(const T * Xi, const T * Xj, size_t n, const double tau, std::vector<double>& ebs, std::vector<std::vector<int>>& weights){
	double eb_Xi = ebs[0];
	double eb_Xj = ebs[1];
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_XiXj = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi, -weights[0][i]), ldexp(eb_Xj, -weights[1][i]));
        double XiXj = Xi[i] * Xj[i];

		error_est_XiXj[i] = e_XiXj;
		error_XiXj[i] = XiXj - XiXj_ori[i];
		if(max_value < error_est_XiXj[i]){
			max_value = error_est_XiXj[i];
			max_index = i;
		}
	}
	// std::cout << "XiXj: max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Xi = ebs[0];
		double eb_Xj = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "coordinate decrease\n";
            double estimate_error_Xi = 0;
            {
                double eb_Xi_ = eb_Xi / 1.5;
                estimate_error_Xi = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi_, -weights[0][i]), ldexp(eb_Xj, -weights[1][i]));
            }
            double estimate_error_Xj = 0;
            {
                double eb_Xj_ = eb_Xj / 1.5;
                estimate_error_Xj = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi, -weights[0][i]), ldexp(eb_Xj_, -weights[1][i]));
            }
			double sum_err = 2 * estimate_error - (estimate_error_Xi + estimate_error_Xj);
            double w_Xi = (estimate_error - estimate_error_Xi) / sum_err;
            double w_Xj = (estimate_error - estimate_error_Xj) / sum_err;
            // === Smooth proportional update ===
            double factor_base = 1.5;
            double alpha = 1.0 - (1.0 / factor_base);
            eb_Xi *= (1.0 - alpha * w_Xi);
            eb_Xj *= (1.0 - alpha * w_Xj);
			estimate_error = compute_bound_multiplication(Xi[i], Xj[i], ldexp(eb_Xi, -weights[0][i]), ldexp(eb_Xj, -weights[1][i]));
		}
		ebs[0] = eb_Xi;
		ebs[1] = eb_Xj;
		return false;
	}
	return true;
}

template<class T>
std::vector<size_t> retrieve_xixj_HPEZ(std::string rdata_file_prefix, std::vector<int> index, T tau, std::vector<T> ebs, size_t num_elements, int weighted, T & max_act_error, T & max_est_error, size_t & weight_file_size, bool decrease_method){
    int max_iter = 30;
	bool tolerance_met = false;
	int n_variable = ebs.size();
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);
    if(!weighted){
        std::vector<PDR::ApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::NegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + species[index[i]];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::HPEZApproximator<T>();
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
            Xi_dec = reconstructed_vars[0].data();
            Xj_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_XiXj = std::vector<T>(num_elements);
            error_est_XiXj = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_XiXj_uniform(Xi_dec, Xj_dec, num_elements, tau, ebs);
            else tolerance_met = halfing_error_XiXj_coordinate(Xi_dec, Xj_dec, num_elements, tau, ebs);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_XiXj);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_XiXj);   	
        }
    }
    else{
        std::vector<PDR::WeightedApproximationBasedReconstructor<T, PDR::HPEZApproximator<T>, MDR::WeightedNegaBinaryBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
        std::vector<std::vector<int>> weights(n_variable, std::vector<int>(num_elements, 0));
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + species[index[i]];
            std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
            std::vector<std::string> files;
            int num_levels = 1;
            for(int i=0; i<num_levels; i++){
                std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
                files.push_back(filename);
            }
            auto approximator = PDR::HPEZApproximator<T>();
            auto encoder = WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto compressor = AdaptiveLevelCompressor(64);
            auto estimator = MaxErrorEstimatorHB<T>();
            auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
            auto retriever = ConcatLevelFileRetriever(metadata_file, files);
            reconstructors.push_back(generateWBPReconstructor<T>(approximator, encoder, compressor, estimator, interpreter, retriever));
            reconstructors.back().fetch_weight = true;
			reconstructors.back().load_metadata();
			weights[i] = reconstructors.back().get_int_weights();
            ebs[i] = ldexp(ebs[i], reconstructors[i].get_max_weight());
        }    
        weight_file_size = reconstructors[0].get_weight_file_size() + reconstructors[1].get_weight_file_size();
        while((!tolerance_met) && (iter < max_iter)){
            iter ++;
            for(int i=0; i<n_variable; i++){
                auto reconstructed_data = reconstructors[i].progressive_reconstruct(ldexp(ebs[i], -reconstructors[i].get_max_weight()), -1);
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
                total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
            }
            Xi_dec = reconstructed_vars[0].data();
            Xj_dec = reconstructed_vars[1].data();
            // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
            // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
            error_XiXj = std::vector<T>(num_elements);
            error_est_XiXj = std::vector<T>(num_elements);
            // std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            if(!decrease_method) tolerance_met = halfing_error_XiXj_uniform(Xi_dec, Xj_dec, num_elements, tau, ebs, weights);
            else tolerance_met = halfing_error_XiXj_coordinate(Xi_dec, Xj_dec, num_elements, tau, ebs, weights);
            // std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
            // MDR::print_vec(ebs);
            /* test
            std::string filename = "./Result/Temp_err.dat";
            std::ofstream outfile1(filename, std::ios::binary);
            if (!outfile1.is_open()) {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                exit(-1);
            }

            outfile1.write(reinterpret_cast<const char*>(error_est_C.data()), error_est_C.size() * sizeof(float));
        
            outfile1.close();
            std::cout << "Data saved successfully to " << filename << std::endl;
            //*/
            // std::cout << names[1] << " requested error = " << tau << std::endl;
            max_act_error = print_max_abs(names[1] + " error", error_XiXj);
            max_est_error = print_max_abs(names[1] + " error_est", error_est_XiXj);   	
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
    int id_i = atoi(argv[argv_id++]);
	int id_j = atoi(argv[argv_id++]);
    if(!(id_i == 1 && id_j == 3) && 
        !(id_i == 4 && id_j == 5) && 
        !(id_i == 0 && id_j == 4) && 
        !(id_i == 3 && id_j == 5)){
        perror("No such QoI\n");
    }

    std::vector<int> index = {id_i, id_j};

    size_t num_elements = 0;
    Xi_ori = MGARD::readfile<T>((data_file_prefix + species[id_i] + ".dat").c_str(), num_elements);
    Xj_ori = MGARD::readfile<T>((data_file_prefix + species[id_j] + ".dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Xi_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Xj_ori)*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Xi_ori, Xj_ori};
    std::vector<double> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    }

	struct timespec start, end;
	int err;
	double elapsed_time;

    std::vector<T> XiXj(num_elements);
	for(int i=0; i<num_elements; i++){
		XiXj[i] = Xi_ori[i] * Xj_ori[i];
	}
	XiXj_ori = XiXj.data();
	double tau = compute_value_range(XiXj)*target_rel_eb;
    T max_act_error = 0, max_est_error = 0;
    size_t weight_file_size = 0;
    std::vector<size_t> total_retrieved_size(n_variable, 0);

    err = clock_gettime(CLOCK_REALTIME, &start);

    switch (compressor)
    {
    case Dummy:
        std::cout << "xixj not applied" << std::endl;
        break;
    case SZ3:
        std::cout << "xixj not applied" << std::endl;
        break;
    case PMGARD:
        std::cout << "xixj not applied" << std::endl;
        break;
    case GE:
        std::cout << "xixj not applied" << std::endl;
        break;
    case HPEZ:
        total_retrieved_size = retrieve_xixj_HPEZ<T>(rdata_file_prefix, index, tau, ebs, num_elements, weighted, max_act_error, max_est_error, weight_file_size, decrease_method);
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