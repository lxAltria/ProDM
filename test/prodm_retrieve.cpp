// prodm_retrieve: the common entry point of ProDM's progressive retrieval pipelines.
//
//   prodm_retrieve <data_dir> <refactor_dir> --vars X,Y,Z --dtype f|d --method mdr|proaicd|pdr|pdr-delta
//                  [--approximator dummy|sz2|sz3|hpez|mgard|ge] [--encoder nega|xor|perbit] [--cp]
//                  [--interpreter greedy|dp|bfs] --tolerance t1 t2 ...
//                  [--qoi <name> --estimator hand --inverse uniform|coordinate] [--joint-range] [--max-iter 30]
//                  [--suffix .dat] [--output <dir>]
//
// Without --qoi, each variable is retrieved to the requested tolerances, given relative
// to its value range, and the error statistics are printed (raw-data error control).
// With --qoi, the variables of the QoI are retrieved together until the QoI error bound
// of the requested tolerance (relative to the QoI value range) is met: the hand-derived
// estimator (SC'24) bounds the QoI error from the per-variable error bounds, and the
// inverse tightens the per-variable bounds (uniformly, or by coordinate descent as in
// HPDC'26) until the estimate meets the tolerance. Weighted refactors (prodm_refactor
// --weights) are detected from the stored weights and handled as in HPDC'26.

#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "ProDM/Decomposer/MultiLevel/MGARDx/utils.hpp"
#include "ProDM/Reconstructor/MDR/Reconstructor.hpp"
#include "ProDM/Reconstructor/PDR/Reconstructor.hpp"
#include "ProDM/Reconstructor/PDR/ApproximationBasedDeltaReconstructor.hpp"
#include "ProDM/Utils/MaskUtils.hpp"

#include "ProDM/Utils/CLI.hpp"
#include "ProDM/QoI/HandQoI.hpp"

using namespace ProDM::CLI;
namespace Hand = ProDM::QoI::Hand;

namespace {

struct Settings {
    Method method;
    Approximator approximator;
    std::string encoder;
    std::string interpreter;
    bool cp;
    std::vector<double> tolerances;
    const Hand::HandQoI* qoi;      // null for raw-data error control
    Hand::Inverse inverse;
    bool coordinate_auto;               // --inverse coordinate: pick the variant from weighted/unweighted
    bool joint_range;
    int max_iter;
    std::string output_dir;
};

void usage(const char* cmd){
    std::cout << "usage: " << cmd << " <data_dir> <refactor_dir> --vars X,Y,Z --dtype f|d --method mdr|proaicd|pdr|pdr-delta\n"
              << "         [--approximator dummy|sz2|sz3|hpez|mgard|ge] [--encoder nega|xor|perbit] [--cp]\n"
              << "         [--interpreter greedy|dp|bfs] --tolerance t1 t2 ...\n"
              << "         [--qoi " << Hand::names() << " --estimator hand --inverse uniform|coordinate|coordinate-min|coordinate-proportional]\n"
              << "         [--joint-range] [--max-iter 30] [--suffix .dat] [--output <dir>]\n"
              << "  tolerances are relative to the value range of each variable (raw-data error control)\n"
              << "  or of the QoI (--qoi); the refactored streams are read from <refactor_dir>/<var>_refactored/\n"
              << "example: " << cmd << " data refactor --vars VelocityX,VelocityY,VelocityZ --dtype d --method pdr --approximator hpez --tolerance 0.01 --qoi Vtot --estimator hand --inverse coordinate" << std::endl;
}

bool file_exists(const std::string& path){
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

// a reconstructor of any pipeline behind a uniform interface
template <class T>
struct Recon {
    std::function<T*(double)> reconstruct;      // progressive reconstruction to an absolute bound
    std::function<size_t()> retrieved_size;
    int max_weight = 0;                         // weighted data: bounds are scaled by 2^max_weight
    size_t weight_file_size = 0;                // stored weights (group-first variable only)
    std::vector<int> int_weights;               // per-point integer weights (weighted data)
};

template <class T, class R>
Recon<T> erase(std::shared_ptr<R> r){
    Recon<T> out;
    out.reconstruct = [r](double tolerance){ return r->progressive_reconstruct(tolerance, -1); };
    out.retrieved_size = [r](){ return r->get_retrieved_size(); };
    return out;
}

template <class T, class F>
void with_approximator(Approximator a, F&& f){
    switch(a){
        case Approximator::DUMMY: f(PDR::DummyApproximator<T>()); return;
#ifdef PRODM_HAVE_SZ2
        case Approximator::SZ2: f(PDR::SZ2Approximator<T>()); return;
#endif
#ifdef PRODM_HAVE_SZ3
        case Approximator::SZ3: f(PDR::SZ3Approximator<T>()); return;
#endif
#ifdef PRODM_HAVE_MGARD
        case Approximator::MGARD: f(PDR::MGARDApproximator<T>()); return;
#endif
#ifdef PRODM_HAVE_HPEZ
        case Approximator::HPEZ: f(PDR::HPEZApproximator<T>()); return;
        case Approximator::GE: f(PDR::GEApproximator<T>()); return;
#endif
        default:
            fail("the requested approximator is not enabled at build time (see the PRODM_WITH_* CMake options)");
    }
}

template <class Estimator, class F>
void with_interpreter(const std::string& name, Estimator estimator, F&& f){
    if(name == "greedy") f(ProDM::SignExcludeGreedyBasedSizeInterpreter<Estimator>(estimator));
    else if(name == "dp") f(ProDM::SignExcludeDPBasedSizeInterpreter<Estimator>(estimator));
    else if(name == "bfs") f(ProDM::SignExcludeBFSBasedSizeInterpreter<Estimator>(estimator));
    else fail("--interpreter must be one of greedy|dp|bfs");
}

// ------------------------------- mdr (SC'21) --------------------------------
template <class T, class T_stream>
Recon<T> make_mdr(const Variable& var, const Settings& s){
    std::string metadata_file = var.refactor_dir + "metadata.bin";
    int num_levels = read_num_levels(metadata_file, false);
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, level_files(var.refactor_dir, num_levels));
    auto estimator = ProDM::MaxErrorEstimatorHB<T>();
    Recon<T> out;
    with_interpreter(s.interpreter, estimator, [&](auto interpreter){
        auto build = [&](auto encoder){
            using R = MDR::ComposedReconstructor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
            auto r = std::make_shared<R>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
            r->load_metadata();
            out = erase<T>(r);
        };
        if(s.encoder == "nega") build(ProDM::NegaBinaryBPEncoder<T, T_stream>());
        else if(s.encoder == "perbit") build(ProDM::PerBitBPEncoder_old<T, T_stream>());
        else fail("--method mdr supports --encoder nega|perbit");
    });
    return out;
}

// ------------------------------ proaicd (SC'26) ------------------------------
template <class T, class T_stream>
Recon<T> make_proaicd(const Variable& var, const Settings& s){
    std::string metadata_file = var.refactor_dir + "metadata.bin";
    std::string data_file = var.refactor_dir + "data.bin";
    int num_dims = 0;
    int num_levels = read_num_levels(metadata_file, false, &num_dims);
    int interpolation_type = read_interpolation_type<T>(metadata_file);
    auto files = level_files(var.refactor_dir, num_levels);
    auto interleaver = MDR::DirectInterleaver_new<T>();
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    Recon<T> out;
    auto dispatch_encoder = [&](auto&& build){
        if(s.encoder == "nega") build(ProDM::NegaBinaryBPEncoder<T, T_stream>());
        else if(s.encoder == "xor") build(ProDM::XORNegaBinaryBPEncoder<T, T_stream>());
        else if(s.encoder == "perbit") build(ProDM::PerBitBPEncoder<T, T_stream>());
        else fail("--method proaicd supports --encoder nega|xor|perbit");
    };
    auto finish = [&](auto r){
        r->load_metadata();
        out = erase<T>(r);
    };
    if(interpolation_type == 0){   // per level
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver<T>();
        auto estimator = MDR::MaxErrorEstimatorHBCubic<T>(num_dims);
        with_interpreter(s.interpreter, estimator, [&](auto interpreter){
            if(s.cp){
                auto retriever = ProDM::OrderedFileRetriever(metadata_file, data_file);
                dispatch_encoder([&](auto encoder){
                    using R = MDR::OrderedCPReconstructor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    finish(std::make_shared<R>(decomposer, interleaver, encoder, compressor, interpreter, retriever));
                });
            } else {
                auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, files);
                dispatch_encoder([&](auto encoder){
                    using R = MDR::FuseComposedReconstructor<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    finish(std::make_shared<R>(decomposer, interleaver, encoder, compressor, interpreter, retriever));
                });
            }
        });
    } else {                        // per layer
        auto decomposer = MDR::MGARDHierarchical_Cubic_Decomposer_Interleaver_new<T>();
        auto estimator = MDR::MaxErrorEstimatorHBCubic_new<T>(1);
        auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, files);
        with_interpreter(s.interpreter, estimator, [&](auto interpreter){
            if(s.cp){
                dispatch_encoder([&](auto encoder){
                    using R = MDR::CPReconstructor_new<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    finish(std::make_shared<R>(decomposer, interleaver, encoder, compressor, interpreter, retriever));
                });
            } else {
                dispatch_encoder([&](auto encoder){
                    using R = MDR::ComposedReconstructor_new<T, decltype(decomposer), decltype(interleaver), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    finish(std::make_shared<R>(decomposer, interleaver, encoder, compressor, interpreter, retriever));
                });
            }
        });
    }
    return out;
}

// --------------------------- pdr (TVCG'23, HPDC'26) --------------------------
// weight bookkeeping shared by the variables of one weight group: the first
// variable of a group loads the stored weights, the others copy them
struct GroupState {
    bool loaded = false;
    std::vector<int> int_weights;
    int max_weight = 0;
};

template <class T>
Recon<T> make_pdr(const Variable& var, const Settings& s, const std::vector<unsigned char>& mask, bool apply_mask, bool weighted, GroupState* gs){
    std::string metadata_file = var.refactor_dir + "metadata.bin";
    auto files = level_files(var.refactor_dir, 1);
    auto compressor = ProDM::AdaptiveLevelCompressor(64);
    auto retriever = ProDM::ConcatLevelFileRetriever(metadata_file, files);
    auto estimator = ProDM::MaxErrorEstimatorHB<T>();
    Recon<T> out;
    with_approximator<T>(s.approximator, [&](auto approximator){
        with_interpreter(s.interpreter, estimator, [&](auto interpreter){
            if(!weighted){
                auto encoder = ProDM::NegaBinaryBPEncoder<T, uint32_t>();
                using R = PDR::ApproximationBasedReconstructor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                auto r = std::make_shared<R>(approximator, encoder, compressor, interpreter, retriever);
                if(apply_mask) r->mask = mask;
                r->load_metadata();
                out = erase<T>(r);
            } else {
                auto encoder = ProDM::WeightedNegaBinaryBPEncoder<T, uint32_t>();
                auto setup = [&](auto r){
                    if(apply_mask) r->mask = mask;
                    if(!gs->loaded){
                        r->fetch_weight = true;
                    } else {
                        r->copy_int_weights(gs->int_weights, gs->max_weight);
                    }
                    r->load_metadata();
                    if(!gs->loaded){
                        gs->int_weights = r->get_int_weights();
                        gs->max_weight = r->get_max_weight();
                        gs->loaded = true;
                        out = erase<T>(r);
                        out.weight_file_size = r->get_weight_file_size();
                    } else {
                        out = erase<T>(r);
                    }
                    out.max_weight = gs->max_weight;
                    out.int_weights = r->get_int_weights();
                };
                if(s.approximator == Approximator::GE){
                    using R = PDR::GEReconstructor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    setup(std::make_shared<R>(approximator, encoder, compressor, interpreter, retriever));
                } else {
                    using R = PDR::WeightedApproximationBasedReconstructor<T, decltype(approximator), decltype(encoder), decltype(compressor), decltype(interpreter), decltype(estimator), decltype(retriever)>;
                    setup(std::make_shared<R>(approximator, encoder, compressor, interpreter, retriever));
                }
            }
        });
    });
    return out;
}

// ------------------------------- pdr-delta ---------------------------------
template <class T>
Recon<T> make_pdr_delta(const Variable& var, const Settings& s){
    std::string dir = var.refactor_dir.substr(0, var.refactor_dir.size() - 1);
    if(!file_exists(dir + "/metadata.bin")) fail("cannot read " + dir + "/metadata.bin; run prodm_refactor first");
    Recon<T> out;
    with_approximator<T>(s.approximator, [&](auto approximator){
        using R = PDR::ApproximationBasedDeltaReconstructor<T, decltype(approximator)>;
        auto r = std::make_shared<R>(approximator, dir);
        r->load_metadata();
        out = erase<T>(r);
    });
    return out;
}

template <class T, class T_stream>
Recon<T> make_reconstructor(const Variable& var, const Settings& s, const std::vector<unsigned char>& mask, bool apply_mask, bool weighted, GroupState* gs){
    if(weighted && s.method != Method::PDR) fail("weighted refactors are supported by --method pdr only");
    switch(s.method){
        case Method::MDR: return make_mdr<T, T_stream>(var, s);
        case Method::PROAICD: return make_proaicd<T, T_stream>(var, s);
        case Method::PDR: return make_pdr<T>(var, s, mask, apply_mask, weighted, gs);
        case Method::PDR_DELTA: return make_pdr_delta<T>(var, s);
    }
    return Recon<T>();
}

// ---------------------------- raw-data error control ----------------------------
template <class T, class T_stream>
void retrieve_raw(const std::vector<Variable>& vars, const Settings& s, const std::string& suffix){
    for(const auto& var : vars){
        if(s.method == Method::PDR && file_exists(var.refactor_dir + "weight.bin")){
            fail(var.name + " was refactored with weights; retrieve it through a QoI (--qoi)");
        }
        std::cout << "Variable " << var.name << std::endl;
        auto data = read_variable<T>(var, 0);
        auto recon = make_reconstructor<T, T_stream>(var, s, {}, false, false, nullptr);
        T value_range = ProDM::compute_value_range(data);
        for(size_t i = 0; i < s.tolerances.size(); i++){
            double tolerance = s.tolerances[i] * value_range;
            std::cout << "Requested tolerance #" << i << " = " << tolerance << " (relative " << s.tolerances[i] << ")" << std::endl;
            double start = now_seconds();
            T* reconstructed = recon.reconstruct(tolerance);
            std::cout << "Reconstruct time: " << now_seconds() - start << "s" << std::endl;
            size_t retrieved_size = recon.retrieved_size();
            std::cout << "Retrieved data size = " << retrieved_size << std::endl;
            ProDM::print_statistics(data.data(), reconstructed, data.size(), retrieved_size);
            std::cout << "Bitrate = " << retrieved_size * 8.0 / data.size() << std::endl;
            if(!s.output_dir.empty()){
                std::string out = s.output_dir + "/" + var.name + suffix;
                ProDM::writefile(out.c_str(), reconstructed, data.size());
            }
            std::cout << std::endl;
        }
    }
}

// ------------------------------ QoI error control ------------------------------
template <class T, class T_stream>
void retrieve_qoi(const std::vector<Variable>& all_vars, const Settings& s, const std::string& refactor_dir){
    using namespace Hand;
    const HandQoI& q = *s.qoi;
    // the variables of the QoI, by canonical index
    std::vector<int> var_index(NUM_VARS, -1);
    for(int v : q.vars){
        var_index[v] = find_variable(all_vars, var_name(v));
        if(var_index[v] < 0) fail(std::string("QoI ") + q.name + " needs variable " + var_name(v) + " (add it to --vars)");
    }
    for(const auto& var : all_vars){
        bool used = false;
        for(int v : q.vars) if(all_vars[var_index[v]].name == var.name) used = true;
        if(!used) std::cout << "note: variable " << var.name << " is not used by QoI " << q.name << " and is not retrieved" << std::endl;
    }

    // originals, value ranges, and exact QoI values
    size_t num_elements = 0;
    std::vector<std::vector<T>> originals(NUM_VARS);
    std::vector<double> ranges(NUM_VARS, 0);
    for(int v : q.vars){
        originals[v] = read_variable<T>(all_vars[var_index[v]], num_elements);
        num_elements = originals[v].size();
        ranges[v] = ProDM::compute_value_range(originals[v]);
    }
    double joint_range = 0;
    {
        T lo = originals[q.vars[0]][0], hi = lo;
        for(int v : q.vars) for(auto x : originals[v]){ lo = std::min(lo, x); hi = std::max(hi, x); }
        joint_range = (double)hi - (double)lo;
    }
    std::vector<double> qoi_ori(num_elements);
    {
        double xi[NUM_VARS] = {0};
        for(size_t i = 0; i < num_elements; i++){
            for(int v : q.vars) xi[v] = originals[v][i];
            qoi_ori[i] = q.value(xi);
        }
    }
    double tau_base = ProDM::compute_value_range(qoi_ori);

    // the velocity mask written by prodm_refactor (SC'24 tools)
    std::vector<unsigned char> mask;
    uint32_t mask_file_size = 0;
    std::string mask_file = refactor_dir + "/mask.bin";
    if(q.uses_velocity_mask && file_exists(mask_file)){
        mask = ProDM::readmask(mask_file.c_str(), mask_file_size);
        if(mask.size() != num_elements) fail("mask " + mask_file + " has " + std::to_string(mask.size()) + " elements, expected " + std::to_string(num_elements));
    }

    // weighted data: the group-first variable holds weight.bin
    bool weighted = false;
    if(s.method == Method::PDR){
        int with = 0, without = 0;
        for(auto g : q.groups){
            int first = group_vars(g)[0];
            if(file_exists(all_vars[var_index[first]].refactor_dir + "weight.bin")) with++; else without++;
        }
        if(with && without) fail("some weight groups of the refactor carry weights and others do not");
        weighted = with > 0;
    }
    Inverse inverse = s.inverse;
    if(s.coordinate_auto) inverse = weighted ? Inverse::COORDINATE_PROPORTIONAL : Inverse::COORDINATE_MIN;
    std::cout << "QoI " << q.name << ", hand estimator, "
              << (inverse == Inverse::UNIFORM ? "uniform" : inverse == Inverse::COORDINATE_MIN ? "coordinate (closest-to-minimum)" : "coordinate (proportional)")
              << " descent, " << (weighted ? "weighted" : "unweighted") << " refactor"
              << (mask.empty() ? "" : ", velocity mask") << std::endl;

    // reconstructors in canonical order (weights are shared inside a group)
    std::vector<Recon<T>> recon(NUM_VARS);
    std::vector<std::vector<int>> weights(NUM_VARS);
    GroupState groups[2];
    size_t weight_file_size = 0;
    for(int v : q.vars){
        GroupState* gs = &groups[v <= VZ ? GROUP_VELOCITY : GROUP_THERMO];
        bool apply_mask = !mask.empty() && v <= VZ;
        recon[v] = make_reconstructor<T, T_stream>(all_vars[var_index[v]], s, mask, apply_mask, weighted, gs);
        weights[v] = recon[v].int_weights;
        weight_file_size += recon[v].weight_file_size;
    }

    double total_start = now_seconds();
    for(double rel : s.tolerances){
        double tau = tau_base * rel;
        std::vector<double> ebs(NUM_VARS, 0);
        for(int v : q.vars){
            ebs[v] = (s.joint_range ? joint_range : ranges[v]) * rel;
            if(weighted) ebs[v] = ldexp(ebs[v], recon[v].max_weight);
        }
        std::vector<const T*> x(NUM_VARS, nullptr);
        std::vector<size_t> retrieved(NUM_VARS, 0);
        int iter = 0;
        bool met = false;
        double max_act_error = 0, max_est_error = 0;
        double start = now_seconds();
        while(!met && iter < s.max_iter){
            iter++;
            for(int v : q.vars){
                double eb = weighted ? ldexp(ebs[v], -recon[v].max_weight) : ebs[v];
                x[v] = recon[v].reconstruct(eb);
                retrieved[v] = recon[v].retrieved_size();
            }
            met = hand_iteration<T>(q, x, num_elements, mask, tau, ebs, weights, inverse, qoi_ori, max_act_error, max_est_error);
            std::cout << "iter " << iter << ": max_est_error = " << max_est_error << ", max_act_error = " << max_act_error << ", next ebs =";
            for(int v : q.vars) std::cout << " " << (weighted ? ldexp(ebs[v], -recon[v].max_weight) : ebs[v]);
            std::cout << std::endl;
        }
        double elapsed = now_seconds() - start;
        size_t total_size = mask_file_size + weight_file_size + std::accumulate(retrieved.begin(), retrieved.end(), size_t(0));
        std::cout << "requested_error = " << tau << std::endl;
        std::cout << "max_est_error = " << max_est_error << std::endl;
        std::cout << "max_act_error = " << max_act_error << std::endl;
        std::cout << "iter = " << iter << std::endl;
        std::cout << "each retrieved size:";
        for(int v : q.vars) std::cout << " " << retrieved[v];
        std::cout << ", mask_file_size = " << mask_file_size << ", weight_file_size = " << weight_file_size << std::endl;
        std::cout << "aggregated cr = " << q.vars.size() * num_elements * sizeof(T) * 1.0 / total_size << std::endl;
        std::cout << "bitrate = " << total_size * 8.0 / (q.vars.size() * num_elements) << std::endl;
        printf("elapsed_time = %.6f\n", elapsed);
        std::cout << std::endl;
    }
    printf("total_elapsed_time = %.6f\n", now_seconds() - total_start);
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
    std::string dtype = opt.choice("dtype", {"f", "d"}, "");
    bool is_double = dtype == "d";
    Settings s;
    s.method = parse_method(opt.get("method"));
    s.approximator = parse_approximator(opt.get("approximator", "hpez"));
    s.encoder = opt.choice("encoder", {"nega", "xor", "perbit"}, "nega");
    s.interpreter = opt.choice("interpreter", {"greedy", "dp", "bfs"}, "greedy");
    s.cp = opt.flag("cp");
    if(opt.flag("no-cp")) s.cp = false;
    for(const auto& t : opt.get_list("tolerance")) s.tolerances.push_back(atof(t.c_str()));
    if(s.tolerances.empty()) fail("--tolerance is required");
    std::string qoi = opt.get("qoi");
    s.qoi = nullptr;
    if(!qoi.empty()){
        s.qoi = Hand::find(qoi);
        if(!s.qoi) fail("unknown QoI '" + qoi + "'; expected one of " + Hand::names());
        if(var_names.empty()) for(int v : s.qoi->vars) var_names.push_back(Hand::var_name(v));
    }
    std::string estimator = opt.choice("estimator", {"hand"}, "hand");
    std::string inverse = opt.choice("inverse", {"uniform", "coordinate", "coordinate-min", "coordinate-proportional"}, "uniform");
    s.coordinate_auto = inverse == "coordinate";
    s.inverse = inverse == "uniform" ? Hand::Inverse::UNIFORM
              : inverse == "coordinate-proportional" ? Hand::Inverse::COORDINATE_PROPORTIONAL
              : Hand::Inverse::COORDINATE_MIN;
    s.joint_range = opt.flag("joint-range");
    s.max_iter = opt.get_int("max-iter", 30);
    s.output_dir = opt.get("output");
    std::string suffix = opt.get("suffix", default_suffix(is_double));
    opt.reject_unknown();
    if(!s.qoi && (opt.has("inverse") || opt.has("joint-range"))) fail("--inverse and --joint-range apply to QoI error control (--qoi)");

    data_dir = absolute_path(data_dir);
    refactor_dir = absolute_path(refactor_dir);
    auto vars = resolve_variables(var_names, data_dir, refactor_dir, suffix);
    if(!s.output_dir.empty()) make_dir(s.output_dir);
    if(is_double){
        if(s.qoi) retrieve_qoi<double, uint64_t>(vars, s, refactor_dir);
        else retrieve_raw<double, uint64_t>(vars, s, suffix);
    } else {
        if(s.qoi) retrieve_qoi<float, uint32_t>(vars, s, refactor_dir);
        else retrieve_raw<float, uint32_t>(vars, s, suffix);
    }
    return 0;
}
