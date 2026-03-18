#ifndef _MDR_BFS_BASED_SIZE_INTERPRETER_HPP
#define _MDR_BFS_BASED_SIZE_INTERPRETER_HPP

#include "SizeInterpreterInterface.hpp"
#include <queue>
#include <unordered_map>
#include <cstdint>
#include <memory>
#include "MDR/RefactorUtils.hpp"

// inorder and round-robin size interpreter

namespace MDR {
    struct State{
        std::vector<uint8_t> pos;
        uint16_t coeff_pos;
        double error;
        uint32_t cost;
        std::shared_ptr<State> parent = nullptr;
    };
    struct PQUnit{
        std::shared_ptr<State> s_ptr;
        bool operator < (const PQUnit& other) const{
            return s_ptr->cost > other.s_ptr->cost;
        }
    };
    // greedy bit-plane retrieval with sign exculsion (excluding the first component)
    template<class ErrorEstimator>
    class SignExcludeBFSBasedSizeInterpreter : public concepts::SizeInterpreterInterface {
    public:
        SignExcludeBFSBasedSizeInterpreter(const ErrorEstimator& e){
            error_estimator = e;
        }
        std::vector<uint32_t> interpret_retrieve_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index) const {
            int num_levels = level_sizes.size();
            int num_bitplanes = level_sizes[0].size();
            std::vector<uint32_t> retrieve_sizes(num_levels, 0);
            // init start state
            current_state->pos = index;
            current_state->error = 0;
            for(int l=0; l<num_levels; l++){
                current_state->error += error_estimator.estimate_error(level_errors[l][index[l]], l, num_levels);
            }
            current_state->cost = 0;
            current_state->parent = nullptr;

            std::priority_queue<PQUnit> pq;
            pq.push({current_state});

            std::vector<std::pair<uint32_t, double>> frontier;
            auto dominated = [&](uint32_t cost, double error){
                for(auto &p : frontier){
                    if(p.first <= cost && p.second <= error)
                        return true;
                }
                return false;
            };

            auto insert_frontier = [&](uint32_t cost, double error){
                std::vector<std::pair<uint32_t,double>> new_frontier;
                for(auto &p : frontier){
                    if(!(cost <= p.first && error <= p.second))
                        new_frontier.push_back(p);
                }
                new_frontier.push_back({cost, error});
                frontier.swap(new_frontier);
            };
            insert_frontier(current_state->cost, current_state->error);

            std::shared_ptr<State> cur;
            while(!pq.empty()){
                cur = pq.top().s_ptr;
                pq.pop();

                // std::cout << "cur->pos:" << std::endl;
                // for(int l=0; l<num_levels; l++){
                //     std::cout << (int)cur->pos[l] << " ";
                // }
                // std::cout << std::endl;
                // std::cout << "cur->error = " << cur->error << ", cur->cost = " << cur->cost << std::endl;
                
                if(cur->error <= tolerance) break;

                uint8_t full_count = 0;
                for(int l=0; l<num_levels; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) full_count++;
                }
                if(full_count == level_sizes.size()) break;

                // Extend states
                for(int l=0; l<num_levels; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) continue;

                    auto nxt = std::make_shared<State>(*cur);

                    nxt->error -= error_estimator.estimate_error(level_errors[l][idx], l, num_levels);
                    nxt->error += error_estimator.estimate_error(level_errors[l][idx + 1], l, num_levels);
                    nxt->cost += level_sizes[l][idx];
                    nxt->parent = cur;
                    nxt->pos[l]++;

                    if(idx != 0 && dominated(nxt->cost, nxt->error)) continue;

                    insert_frontier(nxt->cost, nxt->error);
                    pq.push({nxt});
                }
            }

            for(int l=0; l<num_levels; l++){
                for(int b=index[l]; b<cur->pos[l]; b++){
                    retrieve_sizes[l] += level_sizes[l][b];
                }
            }

            std::cout << "Requested tolerance = " << tolerance << ", estimated error = " << cur->error << std::endl;
            std::cout << "level indexes: " << std::endl;
            for(int i=0; i<num_levels; i++){
                std::cout << (int)cur->pos[i] << " ";
            }
            std::cout << std::endl;
            index = cur->pos;
            _accumulated_error = cur->error;

            std::vector<uint8_t> tmp_path;
            std::vector<double> tmp_error_perstep;

            auto ptr = cur;
            while(ptr->parent != NULL){
                for(int i=0; i<num_levels; i++){
                    if(ptr->parent->pos[i] != ptr->pos[i]) {
                        tmp_path.push_back(static_cast<uint8_t>(i));
                        tmp_error_perstep.push_back(ptr->error);
                    }
                }
                ptr = ptr->parent;
            }
            path.insert(path.end(), tmp_path.rbegin(), tmp_path.rend());
            error_perstep.insert(error_perstep.end(), tmp_error_perstep.rbegin(), tmp_error_perstep.rend());
            // std::reverse(tmp_path.begin(), tmp_path.end());
            // std::reverse(tmp_error_perstep.begin(), tmp_error_perstep.end());
            // // std::cout << "path" << std::endl;
            // for(int i=0; i<tmp_path.size(); i++){
            //     path.push_back(tmp_path[i]);
            //     error_perstep.push_back(tmp_error_perstep[i]);
            //     // std::cout << "(" << (int) tmp_path[i] << ", " << tmp_error_perstep[i] << ") ";
            // }
            // // std::cout << std::endl;

            return retrieve_sizes;
        }

        // For Coeff Prediction Ordering
        std::vector<uint32_t> interpret_coeff_retrieve_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, const uint8_t start_level, const uint8_t num_levels, double tolerance, std::vector<uint8_t>& index) const {
            int end_level = start_level + num_levels;
            std::vector<uint32_t> retrieve_sizes(num_levels, 0);
            // init start state
            current_state->pos = index;
            current_state->error = 0;
            // std::cout << "level_sizes.size() = " << level_sizes.size() << ", level_errors.size() = " << level_errors.size() << std::endl;
            // std::cout << "start_level = " << (int) start_level << ", end_level = " << (int) end_level << std::endl;
            for(int l=start_level; l<end_level; l++){
                current_state->error += error_estimator.estimate_error(level_errors[l][index[l]], l, num_levels);
            }
            // std::cout << "current error = " << current_state->error << std::endl;
            current_state->cost = 0;
            current_state->parent = nullptr;

            std::priority_queue<PQUnit> pq;
            pq.push({current_state});

            std::vector<std::pair<uint32_t, double>> frontier;
            auto dominated = [&](uint32_t cost, double error){
                for(auto &p : frontier){
                    if(p.first <= cost && p.second <= error)
                        return true;
                }
                return false;
            };

            auto insert_frontier = [&](uint32_t cost, double error){
                std::vector<std::pair<uint32_t,double>> new_frontier;
                for(auto &p : frontier){
                    if(!(cost <= p.first && error <= p.second))
                        new_frontier.push_back(p);
                }
                new_frontier.push_back({cost, error});
                frontier.swap(new_frontier);
            };
            insert_frontier(current_state->cost, current_state->error);

            std::shared_ptr<State> cur;
            while(!pq.empty()){
                cur = pq.top().s_ptr;
                pq.pop();

                // std::cout << "cur->pos:" << std::endl;
                // for(int l=0; l<num_levels; l++){
                //     std::cout << (int)cur->pos[l] << " ";
                // }
                // std::cout << std::endl;
                // std::cout << "cur->error = " << cur->error << ", cur->cost = " << cur->cost << std::endl;
                
                if(cur->error <= tolerance) break;

                uint8_t full_count = 0;
                for(int l=start_level; l<end_level; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) full_count++;
                }
                if(full_count == num_levels) break;

                // Extend states
                for(int l=start_level; l<end_level; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) continue;

                    auto nxt = std::make_shared<State>(*cur);

                    nxt->error -= error_estimator.estimate_error(level_errors[l][idx], l, num_levels);
                    nxt->error += error_estimator.estimate_error(level_errors[l][idx + 1], l, num_levels);
                    nxt->cost += level_sizes[l][idx];
                    nxt->parent = cur;
                    nxt->pos[l]++;

                    if(idx != 0 && dominated(nxt->cost, nxt->error)) continue;

                    insert_frontier(nxt->cost, nxt->error);
                    pq.push({nxt});
                }
            }

            // for(int l=start_level; l<end_level; l++){
            //     for(int b=index[l]; b<cur->pos[l]; b++){
            //         retrieve_sizes[l] += level_sizes[l][b];
            //     }
            // }

            // std::cout << "Requested tolerance = " << tolerance << ", estimated error = " << cur->error << std::endl;
            // std::cout << "level indexes: " << std::endl;
            // for(int i=start_level; i<end_level; i++){
            //     std::cout << (int)cur->pos[i] << " ";
            // }
            // std::cout << std::endl;
            index = cur->pos;
            _accumulated_error = cur->error;

            std::vector<uint8_t> tmp_path;
            std::vector<double> tmp_error_perstep;

            auto ptr = cur;
            while(ptr->parent != NULL){
                for(int i=start_level; i<end_level; i++){
                    if(ptr->parent->pos[i] != ptr->pos[i]) {
                        tmp_path.push_back(static_cast<uint8_t>(i));
                        tmp_error_perstep.push_back(ptr->error);
                    }
                }
                ptr = ptr->parent;
            }
            path.insert(path.end(), tmp_path.rbegin(), tmp_path.rend());
            error_perstep.insert(error_perstep.end(), tmp_error_perstep.rbegin(), tmp_error_perstep.rend());
            // std::reverse(tmp_path.begin(), tmp_path.end());
            // std::reverse(tmp_error_perstep.begin(), tmp_error_perstep.end());
            // // std::cout << "path" << std::endl;
            // for(int i=0; i<tmp_path.size(); i++){
            //     path.push_back(tmp_path[i]);
            //     error_perstep.push_back(tmp_error_perstep[i]);
            //     // std::cout << "(" << (int) tmp_path[i] << ", " << tmp_error_perstep[i] << ") ";
            // }
            // // std::cout << std::endl;

            return retrieve_sizes;
        }

        // For Overall Ordering after Coff Prediction Ordering and max ordering
        std::vector<uint32_t> interpret_overall_retrieve_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, const std::vector<uint32_t>& coeff_sizes,
                                                              const std::vector<double>& coeff_error_perstep, const uint8_t coeff_start_level, double tolerance, std::vector<uint8_t>& index, 
                                                              uint16_t & coeff_index) const {
            int num_levels = coeff_start_level + 1;
            std::vector<uint32_t> retrieve_sizes(num_levels, 0);
            // init start state
            current_state->pos = index;
            current_state->coeff_pos = coeff_index;
            current_state->error = 0;
            for(int l=0; l<coeff_start_level; l++){
                current_state->error += error_estimator.estimate_error(level_errors[l][index[l]], l, num_levels);
            }
            current_state->error += coeff_error_perstep[coeff_index];
            current_state->cost = 0;
            current_state->parent = nullptr;
            // std::cout << "input error = " << current_state->error << std::endl; 

            std::priority_queue<PQUnit> pq;
            pq.push({current_state});

            std::vector<std::pair<uint32_t, double>> frontier;
            auto dominated = [&](uint32_t cost, double error){
                for(auto &p : frontier){
                    if(p.first <= cost && p.second <= error)
                        return true;
                }
                return false;
            };

            auto insert_frontier = [&](uint32_t cost, double error){
                std::vector<std::pair<uint32_t,double>> new_frontier;
                for(auto &p : frontier){
                    if(!(cost <= p.first && error <= p.second))
                        new_frontier.push_back(p);
                }
                new_frontier.push_back({cost, error});
                frontier.swap(new_frontier);
            };
            insert_frontier(current_state->cost, current_state->error);

            std::shared_ptr<State> cur;
            while(!pq.empty()){
                cur = pq.top().s_ptr;
                pq.pop();

                // std::cout << "cur->pos:" << std::endl;
                // for(int l=0; l<num_levels; l++){
                //     std::cout << (int)cur->pos[l] << " ";
                // }
                // std::cout << std::endl;
                // std::cout << "cur->error = " << cur->error << ", cur->cost = " << cur->cost << std::endl;
                
                if(cur->error <= tolerance) break;

                uint8_t full_count = 0;
                for(int l=0; l<coeff_start_level; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) full_count++;
                }
                if(cur->coeff_pos >= coeff_sizes.size()) full_count++;
                if(full_count == num_levels) break;

                // Extend states
                for(int l=0; l<coeff_start_level; l++){
                    uint8_t idx = cur->pos[l];
                    if(idx >= level_sizes[l].size()) continue;

                    auto nxt = std::make_shared<State>(*cur);

                    nxt->error -= error_estimator.estimate_error(level_errors[l][idx], l, num_levels);
                    nxt->error += error_estimator.estimate_error(level_errors[l][idx + 1], l, num_levels);
                    nxt->cost += level_sizes[l][idx];
                    nxt->parent = cur;
                    nxt->pos[l]++;

                    if(idx != 0 && dominated(nxt->cost, nxt->error)) continue;

                    insert_frontier(nxt->cost, nxt->error);
                    pq.push({nxt});
                }
                // Extend CP states
                {
                    uint16_t cp_idx = cur->coeff_pos;
                    if(cp_idx >= coeff_sizes.size()) continue;

                    auto nxt = std::make_shared<State>(*cur);

                    nxt->error -= coeff_error_perstep[cp_idx];
                    nxt->error += coeff_error_perstep[cp_idx + 1];
                    nxt->cost += coeff_sizes[cp_idx];
                    nxt->parent = cur;
                    nxt->coeff_pos++;

                    if(cp_idx != 0 && dominated(nxt->cost, nxt->error)) continue;
                    
                    insert_frontier(nxt->cost, nxt->error);
                    pq.push({nxt});
                }
            }

            for(int l=0; l<coeff_start_level; l++){
                for(int b=index[l]; b<cur->pos[l]; b++){
                    retrieve_sizes[l] += level_sizes[l][b];
                }
            }

            std::cout << "Requested tolerance = " << tolerance << ", estimated error = " << cur->error << std::endl;
            std::cout << "level indexes: " << std::endl;
            for(int i=0; i<coeff_start_level; i++){
                std::cout << (int)cur->pos[i] << " ";
            }
            std::cout << (int) cur->coeff_pos << std::endl;
            index = cur->pos;
            coeff_index = cur->coeff_pos;
            _accumulated_error = cur->error;

            std::vector<uint8_t> tmp_path;
            std::vector<double> tmp_error_perstep;

            auto ptr = cur;
            while(ptr->parent != NULL){
                for(int i=0; i<coeff_start_level; i++){
                    if(ptr->parent->pos[i] != ptr->pos[i]) tmp_path.push_back(static_cast<uint8_t>(i));
                }
                if(ptr->parent->coeff_pos != ptr->coeff_pos) tmp_path.push_back(static_cast<uint8_t>(coeff_start_level));
                tmp_error_perstep.push_back(ptr->error);
                ptr = ptr->parent;
            }
            path.insert(path.end(), tmp_path.rbegin(), tmp_path.rend());
            error_perstep.insert(error_perstep.end(), tmp_error_perstep.rbegin(), tmp_error_perstep.rend());
            // std::reverse(tmp_path.begin(), tmp_path.end());
            // std::reverse(tmp_error_perstep.begin(), tmp_error_perstep.end());
            // // std::cout << "path" << std::endl;
            // for(int i=0; i<tmp_path.size(); i++){
            //     path.push_back(tmp_path[i]);
            //     error_perstep.push_back(tmp_error_perstep[i]);
            //     // std::cout << "(" << (int) tmp_path[i] << ", " << tmp_error_perstep[i] << ") ";
            // }
            // // std::cout << std::endl;

            return retrieve_sizes;
        }

        void print() const {
            std::cout << "BFS based size interpreter." << std::endl;
        }
        double get_current_eb(){
            return _accumulated_error;
        }
        void print_path_and_error_perstep(){
            std::cout << path.size() << " " << error_perstep.size() << std::endl;
            for(int i=0; i<path.size(); i++){
                std::cout << "(" << (int)path[i] << ", " << error_perstep[i] << ") ";
            }
            std::cout << std::endl;
        }

        std::vector<uint8_t> get_path(){
            return path;
        }

        std::vector<double> get_error_perstep(){
            return error_perstep;
        }
        void copy_in_cp_levels(std::vector<int> cp_levels){
            error_estimator.cp_levels = cp_levels;
        }
        void reset(){
            path.clear();
            error_perstep.clear();
            current_state = std::make_shared<State>();
            _accumulated_error = 0;
        }
    private:
        ErrorEstimator error_estimator;
        mutable double _accumulated_error = 0;
        mutable std::shared_ptr<State> current_state = std::make_shared<State>();
        mutable std::vector<uint8_t> path;
        mutable std::vector<double> error_perstep;
    };
}
#endif
