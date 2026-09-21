// prodm_refactor: the common entry point of ProDM's refactoring pipelines.
//
//   prodm_refactor <data_dir> <refactor_dir> --vars X,Y,Z --dims n1 n2 n3 --dtype f|d
//                  --method mdr|proaicd|pdr|pdr-delta [--approximator sz3|hpez|ge]
//                  [--target-level 4] [--bitplanes 60] [--eb 1e-3] [--encoder nega|xor|perbit]
//                  [--weights none|hand:<qoi>] [--max-weight 4[,4]] [--block-size 1]
//                  [--prior eb|psnr] [--cp] [--interp-dirs d1 d2 d3] [--mask auto|none] [--suffix .dat]
//
// Each variable <var> is read from <data_dir>/<var><suffix> and refactored into
// <refactor_dir>/<var>_refactored/. The methods route to the pipelines of the
// published papers: mdr (multilevel decomposition, SC'21), proaicd (adaptive
// interpolation and coefficient decomposition, SC'26), pdr (approximation-based
// refactoring, TVCG'23 / HPDC'26; --weights adds the QoI-weighted bitplanes of
// HPDC'26) and pdr-delta (residual snapshots on an approximator).

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Refactor/MDR/Refactor.hpp"
#include "ProDM/Decomposer/MultiLevel/Tuner/Tuner.hpp"
#include "ProDM/Refactor/PDR/Refactor.hpp"
#include "ProDM/Refactor/PDR/ApproximationBasedDeltaRefactor.hpp"
#include "ProDM/Utils/MaskUtils.hpp"

#include "ProDM/Utils/CLI.hpp"
#include "ProDM/QoI/HandQoI.hpp"

using namespace ProDM::CLI;
namespace Hand = ProDM::QoI::Hand;

namespace {

struct Settings {
    Method method;
    Approximator approximator;
    int target_level;
    int num_bitplanes;
    double approximator_eb;
    std::string encoder;
    bool weighted;
    const Hand::HandQoI* weight_qoi;   // recipe of --weights hand:<qoi>
    bool vtot2_recipe;
    std::vector<int> max_weights;           // per active group
    int block_size;
    bool psnr_prior;
    bool cp;
    std::vector<int> interp_dirs;
    bool auto_mask;
};

void usage(const char* cmd){
    std::cout << "usage: " << cmd << " <data_dir> <refactor_dir> --vars X,Y,Z --dims n1 n2 n3 --dtype f|d\n"
              << "         --method mdr|proaicd|pdr|pdr-delta [--approximator sz3|hpez|ge]\n"
              << "         [--target-level 4] [--bitplanes 60] [--eb 1e-3] [--encoder nega|xor|perbit]\n"
              << "         [--weights none|hand:<qoi>] [--max-weight 4[,4]] [--block-size 1]\n"
              << "         [--prior eb|psnr] [--cp] [--interp-dirs d1 d2 d3] [--mask auto|none] [--suffix .dat]\n"
              << "  variables are read from <data_dir>/<var><suffix> (default suffix .dat for d, .dat.f32 for f)\n"
              << "  and written to <refactor_dir>/<var>_refactored/\n"
              << "  --weights hand:<qoi> (pdr only) refactors with the QoI-weighted bitplanes of HPDC'26; qoi is one of "
              << Hand::names() << " or GE (all GE variables)\n"
              << "example: " << cmd << " data refactor --vars VelocityX,VelocityY,VelocityZ --dims 100 500 500 --dtype d --method pdr --approximator hpez" << std::endl;
}

template <class T, class F>
void with_approximator(Approximator a, F&& f){
    switch(a){
#ifdef PRODM_HAVE_SZ3
        case Approximator::SZ3: f(PDR::SZ3Approximator<T>()); return;
#endif
#ifdef PRODM_HAVE_HPEZ
        case Approximator::HPEZ: f(PDR::HPEZApproximator<T>()); return;
        case Approximator::GE: f(PDR::GEApproximator<T>()); return;
#endif
        default:
            fail("the requested approximator is not enabled at build time (see the PRODM_WITH_* CMake options); the command line tools offer sz3, hpez and ge");
    }
}

template <class T, class Refactor>
void timed_refactor(Refactor& refactor, const std::vector<T>& data, const std::vector<uint32_t>& dims, int target_level, int num_bitplanes){
    double start = now_seconds();
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    std::cout << "Refactor time: " << now_seconds() - start << "s" << std::endl;
}

// ------------------------------- mdr (SC'21) --------------------------------
template <class T, class T_stream>
void refactor_mdr(const Variable& var, const std::vector<T>& data, const std::vector<uint32_t>& dims, const Settings& s){
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    auto collector = ProDM::SquaredErrorCollector<T>();
    auto writer = ProDM::ConcatLevelFileWriter(var.refactor_dir + "metadata.bin", level_files(var.refactor_dir, s.target_level + 1));
    auto run = [&](auto encoder, bool negabinary){
        auto refactor = MDR::ComposedRefactor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(collector), decltype(writer)>(decomposer, interleaver, encoder, compressor, collector, writer);
        refactor.negabinary = negabinary;
        timed_refactor(refactor, data, dims, s.target_level, s.num_bitplanes);
    };
    if(s.encoder == "nega") run(ProDM::NegaBinaryBPEncoder<T, T_stream>(), true);
    else if(s.encoder == "perbit") run(ProDM::PerBitBPEncoder_old<T, T_stream>(), false);
    else fail("--method mdr supports --encoder nega|perbit");
}

// ------------------------------ proaicd (SC'26) ------------------------------
template <class T, class T_stream>
void refactor_proaicd(const Variable& var, const std::vector<T>& data, const std::vector<uint32_t>& dims, const Settings& s){
    const int num_dims = dims.size();
    std::string metadata_file = var.refactor_dir + "metadata.bin";
    std::string data_file = var.refactor_dir + "data.bin";
    auto files = level_files(var.refactor_dir, (s.target_level + 2) * num_dims + 1);

    // tune the interpolation order: all-zero means per-level, otherwise per-layer
    std::vector<uint32_t> interp_order;
    {
        auto level_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        auto layer_decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
        auto encoder = ProDM::NegaBinaryBPEncoder<T, T_stream>();
        auto compressor = ProDM::AdaptiveLevelCompressor(64);
        auto level_estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
        auto level_interpreter = ProDM::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic<T>>(level_estimator);
        auto layer_estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
        auto layer_interpreter = ProDM::SignExcludeDPBasedSizeInterpreter<MDR::MaxErrorEstimatorHBCubic_new<T>>(layer_estimator);
        auto tuner = MDR::ProfilingSamplingTuner<T, decltype(level_decomposer), decltype(layer_decomposer), decltype(encoder), decltype(compressor), decltype(level_interpreter), decltype(layer_interpreter), decltype(level_estimator), decltype(layer_estimator)>(level_decomposer, layer_decomposer, encoder, compressor, level_interpreter, layer_interpreter);
        tuner.negabinary = true;
        tuner.mode = s.psnr_prior ? 1 : 0;
        double start = now_seconds();
        tuner.tune(data.data(), dims, s.target_level, s.num_bitplanes, 7, 7);
        std::cout << "Tuner time: " << now_seconds() - start << "s" << std::endl;
        interp_order = tuner.get_best_interp_order();
        std::cout << "best interp order =";
        for(auto o : interp_order) std::cout << " " << o;
        std::cout << std::endl;
    }
    bool per_level = std::accumulate(interp_order.begin(), interp_order.end(), 0u) == 0;

    auto interleaver = MDR::DirectInterleaver_new<T>();
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    auto collector = ProDM::SquaredErrorCollector<T>();
    auto dispatch_encoder = [&](auto&& run){
        if(s.encoder == "nega") run(ProDM::NegaBinaryBPEncoder<T, T_stream>(), true);
        else if(s.encoder == "xor") run(ProDM::XORNegaBinaryBPEncoder<T, T_stream>(), true);
        else if(s.encoder == "perbit") run(ProDM::PerBitBPEncoder<T, T_stream>(), false);
        else fail("--method proaicd supports --encoder nega|xor|perbit");
    };
    if(per_level){
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        if(s.cp){
            auto writer = ProDM::OrderedFileWriter(metadata_file, data_file);
            dispatch_encoder([&](auto encoder, bool negabinary){
                auto refactor = MDR::OrderedCPRefactor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(collector), decltype(writer)>(decomposer, interleaver, encoder, compressor, collector, writer);
                refactor.negabinary = negabinary;
                refactor.coeff_interp_directions = s.interp_dirs;
                refactor.greedy_bfs = true;
                timed_refactor(refactor, data, dims, s.target_level, s.num_bitplanes);
            });
        } else {
            auto writer = ProDM::ConcatLevelFileWriter(metadata_file, files);
            dispatch_encoder([&](auto encoder, bool negabinary){
                auto refactor = MDR::FuseComposedRefactor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(collector), decltype(writer)>(decomposer, interleaver, encoder, compressor, collector, writer);
                refactor.negabinary = negabinary;
                timed_refactor(refactor, data, dims, s.target_level, s.num_bitplanes);
            });
        }
    } else {
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
        decomposer.interp_order = interp_order;
        auto writer = ProDM::ConcatLevelFileWriter(metadata_file, files);
        if(s.cp){
            dispatch_encoder([&](auto encoder, bool negabinary){
                auto refactor = MDR::CPRefactor_new<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(collector), decltype(writer)>(decomposer, interleaver, encoder, compressor, collector, writer);
                refactor.negabinary = negabinary;
                refactor.coeff_interp_directions = s.interp_dirs;
                timed_refactor(refactor, data, dims, s.target_level, s.num_bitplanes);
            });
        } else {
            dispatch_encoder([&](auto encoder, bool negabinary){
                auto refactor = MDR::ComposedRefactor_new<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(collector), decltype(writer)>(decomposer, interleaver, encoder, compressor, collector, writer);
                refactor.negabinary = negabinary;
                timed_refactor(refactor, data, dims, s.target_level, s.num_bitplanes);
            });
        }
    }
}

// --------------------------- pdr (TVCG'23, HPDC'26) --------------------------
// weight bookkeeping shared by the variables of one weight group: the first
// variable of a group derives and stores the integer weights, the others copy them
struct GroupState {
    bool active = false;
    std::vector<int> int_weights;
    bool stored = false;
    int max_weight = 4;
};

template <class T>
void refactor_pdr(const Variable& var, const std::vector<T>& data, const std::vector<uint32_t>& dims, const Settings& s,
                  const std::vector<unsigned char>& mask, bool apply_mask, const std::vector<T>& group_weights, GroupState* gs){
    std::string metadata_file = var.refactor_dir + "metadata.bin";
    auto files = level_files(var.refactor_dir, 1);
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    auto writer = ProDM::ConcatLevelFileWriter(metadata_file, files);
    const uint8_t target_level = 0;
    with_approximator<T>(s.approximator, [&](auto approximator){
        double start = now_seconds();
        if(!s.weighted){
            auto encoder = ProDM::NegaBinaryBPEncoder<T, uint32_t>();
            auto refactor = PDR::ApproximationBasedRefactor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(writer)>(approximator, encoder, compressor, writer);
            refactor.negabinary = true;
            if(apply_mask) refactor.mask = mask;
            refactor.refactor(data.data(), dims, target_level, s.num_bitplanes, (T)s.approximator_eb);
        } else {
            auto encoder = ProDM::WeightedNegaBinaryBPEncoder<T, uint32_t>();
            auto setup = [&](auto& refactor){
                refactor.negabinary = true;
                if(apply_mask) refactor.mask = mask;
                if(!gs->stored){
                    refactor.QoI = group_weights;
                    refactor.store_weight = true;
                } else {
                    refactor.copy_int_weights(gs->int_weights);
                }
            };
            auto finish = [&](auto& refactor){
                if(!gs->stored){
                    gs->int_weights = refactor.get_int_weights();
                    gs->stored = true;
                }
            };
#ifdef PRODM_HAVE_HPEZ
            if constexpr (std::is_same<decltype(approximator), PDR::GEApproximator<T>>::value){
                // the GE variant derives block weights from the mesh block sizes (block_sizes.dat)
                auto refactor = PDR::GERefactor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(writer)>(approximator, encoder, compressor, writer);
                setup(refactor);
                refactor.refactor(data.data(), dims, target_level, s.num_bitplanes, (T)s.approximator_eb, gs->max_weight);
                finish(refactor);
            } else
#endif
            {
                auto refactor = PDR::WeightedApproximationBasedRefactor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(writer)>(approximator, encoder, compressor, writer);
                setup(refactor);
                refactor.refactor(data.data(), dims, target_level, s.num_bitplanes, (T)s.approximator_eb, gs->max_weight, s.block_size);
                finish(refactor);
            }
        }
        std::cout << "Refactor time: " << now_seconds() - start << "s" << std::endl;
    });
}

// ------------------------------- pdr-delta ---------------------------------
template <class T>
void refactor_pdr_delta(const Variable& var, const std::vector<T>& data, const std::vector<uint32_t>& dims, const Settings& s){
    std::string dir = var.refactor_dir.substr(0, var.refactor_dir.size() - 1);
    with_approximator<T>(s.approximator, [&](auto approximator){
        auto refactor = PDR::ApproximationBasedDeltaRefactor<T, decltype(approximator)>(approximator, dir);
        // target level and bitplanes are not used by the delta refactor
        timed_refactor(refactor, data, dims, 0, 0);
    });
}

template <class T>
void run(const std::vector<Variable>& vars, const std::vector<uint32_t>& dims, const std::string& refactor_dir, Settings s){
    using T_stream = typename std::conditional<std::is_same<T, float>::value, uint32_t, uint64_t>::type;
    const int max_bitplanes = std::is_same<T, float>::value ? 32 : (s.method == Method::PDR ? 60 : 64);
    if(s.num_bitplanes <= 0) s.num_bitplanes = max_bitplanes == 64 ? 60 : max_bitplanes;
    if(s.num_bitplanes % 2 == 1){
        s.num_bitplanes++;
        std::cout << "Change to " << s.num_bitplanes << " bitplanes for simplicity of negabinary encoding" << std::endl;
    }
    if(s.num_bitplanes > max_bitplanes){
        s.num_bitplanes = max_bitplanes;
        std::cout << "At most " << max_bitplanes << " bitplanes are supported for this data type and method" << std::endl;
    }
    if(s.method == Method::PROAICD){
        int min_dim = *std::min_element(dims.begin(), dims.end()) - 1;
        int max_level = log2(min_dim);
        if(s.target_level > max_level){
            s.target_level = max_level;
            std::cout << "Target level exceeds the min dimension, change to target_level = " << s.target_level << std::endl;
        }
    }

    size_t num_elements = 1;
    for(auto d : dims) num_elements *= d;
    std::vector<std::vector<T>> data(vars.size());
    for(size_t i = 0; i < vars.size(); i++) data[i] = read_variable<T>(vars[i], num_elements);

    // the velocity mask (nonzero total velocity) when the three velocities are present
    std::vector<unsigned char> mask;
    std::vector<int> vel_index(3, -1);
    bool have_velocities = true;
    for(int k = 0; k < 3; k++){
        vel_index[k] = find_variable(vars, velocity_names()[k]);
        if(vel_index[k] < 0) have_velocities = false;
    }
    // the mask is applied by the pdr encoders only; the other pipelines encode every point
    if(s.auto_mask && have_velocities && s.method == Method::PDR){
        mask.assign(num_elements, 0);
        size_t valid = 0;
        for(size_t i = 0; i < num_elements; i++){
            T vx = data[vel_index[0]][i], vy = data[vel_index[1]][i], vz = data[vel_index[2]][i];
            if(vx*vx + vy*vy + vz*vz != 0){ mask[i] = 1; valid++; }
        }
        std::string mask_file = refactor_dir + "/mask.bin";
        ProDM::writemask(mask_file.c_str(), mask.data(), mask.size());
        std::cout << "velocity mask: " << valid << " of " << num_elements << " points nonzero, written to " << mask_file << std::endl;
    }

    // weight groups of --weights hand:<qoi>
    GroupState groups[2];
    std::vector<int> var_group(vars.size(), -1);
    std::vector<std::vector<T>> group_weights(2);
    if(s.weighted){
        if(s.method != Method::PDR) fail("--weights is supported by --method pdr only");
        std::vector<const T*> x(Hand::NUM_VARS, nullptr);
        for(int v = 0; v < Hand::NUM_VARS; v++){
            int idx = find_variable(vars, Hand::var_name(v));
            if(idx >= 0) x[v] = data[idx].data();
        }
        int active = 0;
        for(auto g : s.weight_qoi->groups){
            bool complete = true;
            for(int v : Hand::group_vars(g)) if(!x[v]) complete = false;
            if(!complete) continue;
            groups[g].active = true;
            groups[g].max_weight = s.max_weights[std::min<size_t>(active, s.max_weights.size() - 1)];
            group_weights[g] = Hand::group_weights<T>(g, s.vtot2_recipe, x, num_elements, mask);
            for(int v : Hand::group_vars(g)) var_group[find_variable(vars, Hand::var_name(v))] = g;
            std::cout << "weight group " << Hand::group_name(g) << ": max weight " << groups[g].max_weight << std::endl;
            active++;
        }
        if(active == 0) fail("--weights hand:" + s.weight_qoi->name + " needs the variables of at least one of its groups (VelocityX/Y/Z or Pressure/Density)");
        for(size_t i = 0; i < vars.size(); i++){
            if(var_group[i] < 0) fail("variable " + vars[i].name + " is not covered by the weights of hand:" + s.weight_qoi->name + "; refactor it separately with --weights none");
        }
    }

    // weighted refactoring stores a group's weights with its canonical-first variable (VelocityX, Pressure),
    // so those are processed first whatever the order of --vars; unweighted runs keep the given order
    std::vector<size_t> order(vars.size());
    for(size_t i = 0; i < vars.size(); i++) order[i] = i;
    if(s.weighted){
        auto rank = [&](size_t i){
            for(int v = 0; v < Hand::NUM_VARS; v++) if(vars[i].name == Hand::var_name(v)) return v;
            return Hand::NUM_VARS + (int)i;
        };
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b){ return rank(a) < rank(b); });
    }

    double start = now_seconds();
    for(size_t i : order){
        make_dir(vars[i].refactor_dir);
        std::cout << "Refactoring " << vars[i].name << " -> " << vars[i].refactor_dir << std::endl;
        bool apply_mask = !mask.empty() && is_velocity(vars[i].name);
        switch(s.method){
            case Method::MDR: refactor_mdr<T, T_stream>(vars[i], data[i], dims, s); break;
            case Method::PROAICD: refactor_proaicd<T, T_stream>(vars[i], data[i], dims, s); break;
            case Method::PDR: {
                int g = var_group[i];
                refactor_pdr<T>(vars[i], data[i], dims, s, mask, apply_mask, g < 0 ? std::vector<T>() : group_weights[g], g < 0 ? nullptr : &groups[g]);
                break;
            }
            case Method::PDR_DELTA: refactor_pdr_delta<T>(vars[i], data[i], dims, s); break;
        }
    }
    printf("elapsed_time = %.6f\n", now_seconds() - start);
}

}

int main(int argc, char** argv){
    if(argc < 3){
        usage(argv[0]);
        return 0;
    }
    Options opt(argc, argv, 2);
    if(opt.positional.size() != 2) fail("expected <data_dir> <refactor_dir> before the options");
    std::string data_dir = opt.positional[0];
    std::string refactor_dir = opt.positional[1];

    auto var_names = opt.get_list("vars");
    auto dim_strs = opt.get_list("dims");
    std::vector<uint32_t> dims;
    for(const auto& d : dim_strs) dims.push_back(atoi(d.c_str()));
    if(dims.empty()) fail("--dims is required");
    for(auto d : dims) if(d == 0) fail("--dims must be positive integers");
    std::string dtype = opt.choice("dtype", {"f", "d"}, "");
    bool is_double = dtype == "d";

    Settings s;
    s.method = parse_method(opt.get("method"));
    s.approximator = parse_approximator(opt.get("approximator", "hpez"));
    s.target_level = opt.get_int("target-level", 4);
    s.num_bitplanes = opt.get_int("bitplanes", 0);
    s.approximator_eb = opt.get_double("eb", 0.001);
    s.encoder = opt.choice("encoder", {"nega", "xor", "perbit"}, "nega");
    s.block_size = opt.get_int("block-size", 1);
    s.psnr_prior = opt.choice("prior", {"eb", "psnr"}, "eb") == "psnr";
    s.cp = opt.flag("cp");
    if(opt.flag("no-cp")) s.cp = false;
    for(const auto& d : opt.get_list("interp-dirs")) s.interp_dirs.push_back(atoi(d.c_str()));
    if(!s.interp_dirs.empty() && s.interp_dirs.size() != dims.size()) fail("--interp-dirs needs one entry per dimension");
    s.auto_mask = opt.choice("mask", {"auto", "none"}, "auto") == "auto";
    std::string suffix = opt.get("suffix", default_suffix(is_double));

    std::string weights = opt.get("weights", "none");
    s.weighted = weights != "none";
    s.weight_qoi = nullptr;
    s.vtot2_recipe = false;
    if(s.weighted){
        if(weights.rfind("hand:", 0) != 0) fail("--weights must be none or hand:<qoi>");
        std::string qoi = weights.substr(5);
        if(qoi == "GE") qoi = "Mach";   // both GE groups
        s.weight_qoi = Hand::find(qoi);
        if(!s.weight_qoi) fail("unknown QoI '" + qoi + "' in --weights; expected one of " + Hand::names() + " or GE");
        s.vtot2_recipe = qoi == "Vtot2";
    }
    auto mw = opt.get_list("max-weight");
    if(mw.empty()) s.max_weights = {4};
    for(const auto& m : mw) s.max_weights.push_back(atoi(m.c_str()));
    opt.reject_unknown();

    if((s.method == Method::PDR || s.method == Method::PDR_DELTA) && !opt.has("approximator")){
        std::cout << "using the default approximator hpez (set --approximator to change)" << std::endl;
    }

    make_dir(refactor_dir);
    data_dir = absolute_path(data_dir);
    refactor_dir = absolute_path(refactor_dir);
    auto vars = resolve_variables(var_names, data_dir, refactor_dir, suffix);
    if(is_double) run<double>(vars, dims, refactor_dir, s);
    else run<float>(vars, dims, refactor_dir, s);
    return 0;
}
