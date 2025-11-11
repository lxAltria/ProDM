#ifndef _MDR_DP_BASED_SIZE_INTERPRETER_HPP
#define _MDR_DP_BASED_SIZE_INTERPRETER_HPP

#include "SizeInterpreterInterface.hpp"
#include <queue>
#include "MDR/RefactorUtils.hpp"

// inorder and round-robin size interpreter

namespace MDR {
    // greedy bit-plane retrieval with sign exculsion (excluding the first component)
    template<class ErrorEstimator>
    class SignExcludeDPBasedSizeInterpreter : public concepts::SizeInterpreterInterface {
    public:
        SignExcludeDPBasedSizeInterpreter(const ErrorEstimator& e){
            error_estimator = e;
        }
        std::vector<uint32_t> interpret_retrieve_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index) const {
            int num_levels = level_sizes.size();
            int num_bitplanes = level_sizes[0].size();
            if(first_time){
                // std::vector<std::vector<double>> err(num_levels, std::vector<double>(num_bitplanes, 0));
                err = std::vector<std::vector<double>>(num_levels, std::vector<double>(num_bitplanes, 0));
                for(int l=0; l<num_levels; l++){
                    for(int b=0; b<num_bitplanes; b++){
                        err[l][b] = error_estimator.estimate_error(level_errors[l][b], l, num_levels);
                    }
                }
                // std::vector<std::vector<uint32_t>> SavedSize(num_levels, std::vector<uint32_t>(num_bitplanes + 1, 0));
                SavedSize = std::vector<std::vector<uint32_t>>(num_levels, std::vector<uint32_t>(num_bitplanes + 1, 0));
                for(int l=0; l<num_levels; l++){ 
                    // SavedSize[l][0] = 0;
                    for(int b=1; b<num_bitplanes+1; b++){
                        int bp = num_bitplanes - b;
                        SavedSize[l][b] = SavedSize[l][b-1] + level_sizes[l][bp];
                    }
                }
                first_time = false;
            }
            // int ERR_MIN = 0;
            int ERR_MAX = 1e3;
            // int ERR_RANGE = ERR_MAX - ERR_MIN;
            std::vector<std::vector<int>> cost(num_levels, std::vector<int>(num_bitplanes + 1));
            for(int l=0; l<num_levels; l++){
                cost[l][0] = 0;
                // std::cout << "level " << l << " costs: " << std::endl;
                for (int b = 1; b < num_bitplanes + 1; ++b) {
                    int bp = num_bitplanes - b;
                    double ratio = std::max(0.0, err[l][bp] / (tolerance / ERR_MAX)); // err[l][bp]/1e-6
                    // double scaled = ratio * ERR_MAX;
                    ratio = std::min(ratio, double(ERR_MAX));
                    int c = (int)std::ceil(ratio);
                    // if (c < ERR_MIN) c = ERR_MIN;
                    // if (c > ERR_MAX) c = ERR_MAX;
                    cost[l][b] = c;
                    // cost[l][b] = err[l][bp];
                    // std::cout << "cost[" << l << "][" << b << "] = " << cost[l][b] << std::endl;
                    // std::cout << cost[l][b] << " ";
                }
                // std::cout << std::endl;
            }

            // const int NEG_INF = -1000000000;
            std::vector<std::vector<int>> DP(num_levels + 1, std::vector<int>(ERR_MAX + 1, 0));
            DP[0][0] = 0;
            std::vector<std::vector<int>> choose(num_levels, std::vector<int>(ERR_MAX + 1, -1));

            for(int l=1; l<num_levels+1; l++){
                for(int e=0; e<ERR_MAX+1; e++){
                    int best_b = 0;
                    int best_val = DP[l-1][e];
                    for(int b=1; b<num_bitplanes+1; b++){
                        int w = cost[l-1][b];
                        int v = SavedSize[l-1][b];
                        if(e >= w){// && DP[l-1][e-w] != NEG_INF){
                            // if((e-w <= 0) || (e-w >= ERR_MAX + 1)) std::cout << "l-1 = " << l-1 << ", e-w = " << e - w << std::endl;
                            int cand = DP[l-1][e-w] + v;
                            if(cand > best_val){
                                best_val = cand;
                                best_b = b;
                            }
                        }
                    }
                    DP[l][e] = best_val;
                    choose[l-1][e] = best_b;
                }
            }

            int best_e = 0;
            for(int e=0; e<ERR_MAX; e++){
                if(DP[num_levels][e] > DP[num_levels][best_e]) best_e = e;
            }
            // std::cout << "best_e = " << best_e << std::endl;
            // uint32_t best_saved = DP[num_levels][best_e];

            std::vector<int> dropped(num_levels, 0);
            int e = best_e;
            for(int l=num_levels - 1; l>=0; l--){
                dropped[l] = choose[l][e];
                // if(dropped[l] > cost[0].size()) std::cout << "dropped[l] = " << dropped[l] << std::endl;
                e -= cost[l][dropped[l]];
            }
            std::vector<uint32_t> retrieve_sizes(num_levels, 0);
            double accumulated_error = 0;
            // std::cout << "level indexes: " << std::endl;
            for(int l=0; l<num_levels; l++){
                uint8_t new_index = uint8_t(num_bitplanes - dropped[l]);  // drop 12, 32-12=20 0~19
                for(int b=index[l]; b < new_index; b++){
                    retrieve_sizes[l] += level_sizes[l][b];
                }
                index[l] = (new_index > index[l]) ? new_index : index[l];
                // std::cout << (int)index[l] << " ";
                if(index[l] < num_bitplanes)
                    accumulated_error += err[l][index[l]];
            }
            // std::cout << std::endl;
            _accumulated_error = accumulated_error;
            // std::cout << "Tolerance = " << tolerance << ", accumulated_error = " << _accumulated_error << std::endl;
            return retrieve_sizes;
        }
        void print() const {
            std::cout << "Greedy based size interpreter." << std::endl;
        }
        double get_current_eb(){
            return _accumulated_error;
        }
    private:
        ErrorEstimator error_estimator;
        mutable double _accumulated_error = 0;
        mutable std::vector<std::vector<double>> err;
        mutable std::vector<std::vector<uint32_t>> SavedSize;
        mutable bool first_time = true;
    };
}
#endif
