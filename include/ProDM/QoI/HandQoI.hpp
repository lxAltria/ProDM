#ifndef PRODM_QOI_HANDQOI_HPP
#define PRODM_QOI_HANDQOI_HPP

// QoI/HandQoI.hpp - the QoI registry of prodm_retrieve / prodm_refactor with the hand-derived error
// estimators of the SC'24 method (Wu et al.) and the error-bound descents of the
// SC'24 / HPDC'26 tools (app/GE/qoi_*_d64.cpp and the former test_qoi_reconstructor).
//
// A QoI is a formula over the canonical variables VelocityX, VelocityY, VelocityZ,
// Pressure, Density. For each QoI the registry gives its value and the hand bound
// built from the primitive bounds in ProDM/Legacy/QoIUtils.hpp; the descents
// (uniform, coordinate) are written once over that bound. The weight recipes of
// the HPDC'26 refactoring are also listed here (1/V_total on the velocities,
// 1/Density^2 on the thermodynamic variables).

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "ProDM/Legacy/QoIUtils.hpp"

#include "ProDM/Namespace.hpp"

namespace ProDM::QoI::Hand {

enum Var { VX = 0, VY = 1, VZ = 2, P = 3, D = 4, NUM_VARS = 5 };

inline const char* var_name(int v){
    static const char* names[NUM_VARS] = {"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};
    return names[v];
}

// gas constants of the GE QoIs
constexpr double R_gas = 287.1;
constexpr double gamma_ = 1.4;
constexpr double mu_r = 1.716e-5;
constexpr double T_r = 273.15;
constexpr double S_c = 110.4;

using ProDM::Legacy::compute_bound_x_square;
using ProDM::Legacy::compute_bound_square_root_x;
using ProDM::Legacy::compute_bound_radical;
using ProDM::Legacy::compute_bound_multiplication;
using ProDM::Legacy::compute_bound_division;

// x[] and eb[] are indexed by canonical variable (NUM_VARS entries; unused ones are ignored).
// masked == false means the point has zero velocity and its velocity terms are exact.
// Each QoI provides value(x) and eval(x, eb, masked, value), which returns the hand
// bound and the value from one pass over the formula (the bound chains are those of
// app/GE/qoi_*_d64.cpp and the former test_qoi_reconstructor).

// --- V_total -----------------------------------------------------------------
struct Vtot {
    static double value(const double* x){
        return sqrt(x[VX]*x[VX] + x[VY]*x[VY] + x[VZ]*x[VZ]);
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        double V2 = x[VX]*x[VX] + x[VY]*x[VY] + x[VZ]*x[VZ];
        value = sqrt(V2);
        if(!masked) return 0;
        double e_V2 = compute_bound_x_square(x[VX], eb[VX]) + compute_bound_x_square(x[VY], eb[VY]) + compute_bound_x_square(x[VZ], eb[VZ]);
        return compute_bound_square_root_x(V2, e_V2);
    }
};

// --- V_total^2 ---------------------------------------------------------------
struct Vtot2 {
    static double value(const double* x){
        return x[VX]*x[VX] + x[VY]*x[VY] + x[VZ]*x[VZ];
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        value = x[VX]*x[VX] + x[VY]*x[VY] + x[VZ]*x[VZ];
        if(!masked) return 0;
        return compute_bound_x_square(x[VX], eb[VX]) + compute_bound_x_square(x[VY], eb[VY]) + compute_bound_x_square(x[VZ], eb[VZ]);
    }
};

// --- temperature T = P / (D R) ---------------------------------------------------
struct Temperature {
    static double value(const double* x){
        return x[P] / (x[D] * R_gas);
    }
    static double eval(const double* x, const double* eb, bool, double& value){
        value = x[P] / (x[D] * R_gas);
        return (1.0 / R_gas) * compute_bound_division(x[P], x[D], eb[P], eb[D]);
    }
};

// --- speed of sound C = sqrt(gamma R T) -----------------------------------------
struct SoundSpeed {
    static double value(const double* x){
        return sqrt(gamma_ * R_gas) * sqrt(Temperature::value(x));
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        const double c_2 = sqrt(gamma_ * R_gas);
        double Temp;
        double e_T = Temperature::eval(x, eb, masked, Temp);
        value = c_2 * sqrt(Temp);
        return c_2 * compute_bound_square_root_x(Temp, e_T);
    }
};

// --- Mach = V_total / C ------------------------------------------------------------
struct Mach {
    static double value(const double* x){
        return Vtot::value(x) / SoundSpeed::value(x);
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        double V, C;
        double e_V = Vtot::eval(x, eb, masked, V);
        double e_C = SoundSpeed::eval(x, eb, masked, C);
        value = V / C;
        return compute_bound_division(V, C, e_V, e_C);
    }
};

// --- total pressure PT = P (1 + (gamma-1)/2 Mach^2)^3.5 ------------------------------
struct TotalPressure {
    static double value(const double* x){
        double M = Mach::value(x);
        double mt = 1 + (gamma_ - 1) / 2 * M * M;
        return x[P] * sqrt(pow(mt, 7));
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        static const int C7[8] = {1, 7, 21, 35, 35, 21, 7, 1};
        double M;
        double e_M = Mach::eval(x, eb, masked, M);
        double e_mt = ldexp(gamma_ - 1, -1) * compute_bound_x_square(M, e_M);
        double mt = 1 + ldexp(gamma_ - 1, -1) * M * M;
        // error of mt^7 by the binomial expansion, fed directly into the final product
        // (the square root is not applied to it, as in the SC'24 implementation)
        double mt_pow[8], e_mt_pow[8];
        mt_pow[0] = 1; e_mt_pow[0] = 1;
        for(int k = 1; k <= 7; k++){
            mt_pow[k] = mt_pow[k-1] * mt;
            e_mt_pow[k] = e_mt_pow[k-1] * e_mt;
        }
        double e_mt_mi = 0;
        for(int k = 1; k <= 7; k++) e_mt_mi += C7[k] * mt_pow[7-k] * e_mt_pow[k];
        double mt_mi = sqrt(mt_pow[7]);
        value = x[P] * mt_mi;
        return compute_bound_multiplication(x[P], mt_mi, eb[P], e_mt_mi);
    }
};

// --- viscosity mu (Sutherland's law) ----------------------------------------------
struct Viscosity {
    static double value(const double* x){
        double Temp = Temperature::value(x);
        double TrS_TS = (T_r + S_c) / (Temp + S_c);
        double T_Tr_3 = pow(Temp / T_r, 3);
        return mu_r * sqrt(T_Tr_3) * TrS_TS;
    }
    static double eval(const double* x, const double* eb, bool masked, double& value){
        const double c_3 = T_r + S_c;
        double Temp;
        double e_T = Temperature::eval(x, eb, masked, Temp);
        double e_TrS_TS = c_3 * compute_bound_radical(Temp, S_c, e_T);
        double TrS_TS = c_3 / (Temp + S_c);
        double e_T_Tr_3 = 3*(Temp/T_r)*(Temp/T_r)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
        double T_Tr_3 = (Temp/T_r)*(Temp/T_r)*(Temp/T_r);
        double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
        double T_Tr_3_sqrt = sqrt(T_Tr_3);
        value = mu_r * T_Tr_3_sqrt * TrS_TS;
        return mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
    }
};

// --- registry --------------------------------------------------------------------
enum Group { GROUP_VELOCITY = 0, GROUP_THERMO = 1 };

enum Id { ID_VTOT, ID_VTOT2, ID_T, ID_C, ID_MACH, ID_PT, ID_MU };

struct HandQoI {
    std::string name;
    Id id;
    std::vector<int> vars;               // canonical variables the QoI reads (also the ones retrieved)
    std::vector<Group> groups;           // weight groups of the HPDC'26 recipe
    bool uses_velocity_mask;             // velocity terms vanish at zero-velocity points
    double (*value)(const double*);      // exact value (used for the original data)
};

inline const std::vector<HandQoI>& registry(){
    static const std::vector<HandQoI> qois = {
        {"Vtot",  ID_VTOT,  {VX, VY, VZ},       {GROUP_VELOCITY},               true,  Vtot::value},
        {"Vtot2", ID_VTOT2, {VX, VY, VZ},       {GROUP_VELOCITY},               true,  Vtot2::value},
        {"T",     ID_T,     {P, D},             {GROUP_THERMO},                 false, Temperature::value},
        {"C",     ID_C,     {P, D},             {GROUP_THERMO},                 false, SoundSpeed::value},
        {"Mach",  ID_MACH,  {VX, VY, VZ, P, D}, {GROUP_VELOCITY, GROUP_THERMO}, true,  Mach::value},
        {"PT",    ID_PT,    {VX, VY, VZ, P, D}, {GROUP_VELOCITY, GROUP_THERMO}, true,  TotalPressure::value},
        {"mu",    ID_MU,    {P, D},             {GROUP_THERMO},                 false, Viscosity::value},
    };
    return qois;
}

// call f<Q>() with the QoI type of q, so that the per-point evaluation inlines
template <class F>
auto visit(const HandQoI& q, F&& f){
    switch(q.id){
        case ID_VTOT:  return f(Vtot{});
        case ID_VTOT2: return f(Vtot2{});
        case ID_T:     return f(Temperature{});
        case ID_C:     return f(SoundSpeed{});
        case ID_MACH:  return f(Mach{});
        case ID_PT:    return f(TotalPressure{});
        case ID_MU:    return f(Viscosity{});
    }
    return f(Vtot{});
}

inline const HandQoI* find(const std::string& name){
    for(const auto& q : registry()) if(q.name == name) return &q;
    return nullptr;
}

inline std::string names(){
    std::string s;
    for(const auto& q : registry()) s += (s.empty() ? "" : "|") + q.name;
    return s;
}

inline const char* group_name(Group g){
    return g == GROUP_VELOCITY ? "velocity" : "thermodynamic";
}

inline std::vector<int> group_vars(Group g){
    return g == GROUP_VELOCITY ? std::vector<int>{VX, VY, VZ} : std::vector<int>{P, D};
}

// per-point weights of a group (HPDC'26): 1/V_total on the velocities (0 where the
// mask is 0), |Vx|+|Vy|+|Vz| for the V_total^2 recipe, 1/Density^2 on Pressure/Density
template <class T>
std::vector<T> group_weights(Group g, bool vtot2_recipe, const std::vector<const T*>& x, size_t n, const std::vector<unsigned char>& mask){
    std::vector<T> w(n, 0);
    if(g == GROUP_VELOCITY){
        for(size_t i = 0; i < n; i++){
            if(vtot2_recipe){
                w[i] = fabs(x[VX][i]) + fabs(x[VY][i]) + fabs(x[VZ][i]);
            } else if(mask.empty() || mask[i]){
                T V = sqrt(x[VX][i]*x[VX][i] + x[VY][i]*x[VY][i] + x[VZ][i]*x[VZ][i]);
                w[i] = 1.0 / V;
            }
        }
    } else {
        for(size_t i = 0; i < n; i++) w[i] = 1.0 / (x[D][i] * x[D][i]);
    }
    return w;
}

// --- error-bound descent ------------------------------------------------------------
enum class Inverse { UNIFORM, COORDINATE_MIN, COORDINATE_PROPORTIONAL };

// One retrieval iteration of the SC'24 loop: evaluate the hand bound on every point
// with the current per-variable bounds ebs[] (weighted points use ldexp(eb, -w)), and,
// if the maximum estimate exceeds tau, tighten ebs[] at the worst point until the
// estimate there meets tau. Returns true when the tolerance is met. x[] holds the
// reconstructed variables by canonical index (unused ones may be null), qoi_ori the
// exact QoI values; weights[] (by canonical index) may be empty for unweighted data.
template <class T, class Q>
bool hand_iteration_impl(Q, const HandQoI& q, const std::vector<const T*>& x, size_t n, const std::vector<unsigned char>& mask,
                         double tau, std::vector<double>& ebs, const std::vector<std::vector<int>>& weights, Inverse inverse,
                         const std::vector<double>& qoi_ori, double& max_act_error, double& max_est_error){
    const bool masked_qoi = q.uses_velocity_mask && !mask.empty();
    const int nv = q.vars.size();
    const int* vars = q.vars.data();
    bool weighted = false;
    for(int k = 0; k < nv; k++) if(!weights[vars[k]].empty()) weighted = true;
    const unsigned char* mask_ptr = masked_qoi ? mask.data() : nullptr;
    const int* wptr[NUM_VARS] = {nullptr};
    const T* xptr[NUM_VARS] = {nullptr};
    for(int k = 0; k < nv; k++){
        xptr[vars[k]] = x[vars[k]];
        if(weighted && !weights[vars[k]].empty()) wptr[vars[k]] = weights[vars[k]].data();
    }
    // effective bounds of the current point: fixed when the data is unweighted; the
    // weighted scaling 2^-w is taken from a table (exact, as ldexp) instead of calling ldexp per point
    std::vector<double> scale;
    if(weighted){
        int max_w = 0;
        for(int k = 0; k < nv; k++){
            if(!wptr[vars[k]]) continue;
            for(size_t i = 0; i < n; i++) if(wptr[vars[k]][i] > max_w) max_w = wptr[vars[k]][i];
        }
        scale.resize(max_w + 1);
        for(int w = 0; w <= max_w; w++) scale[w] = ldexp(1.0, -w);
    }
    auto effective_eb = [&](size_t i, const std::vector<double>& e, double* ei){
        for(int k = 0; k < nv; k++){
            int v = vars[k];
            ei[v] = wptr[v] ? e[v] * scale[wptr[v][i]] : e[v];
        }
    };
    double xi[NUM_VARS] = {0}, ei[NUM_VARS] = {0};
    if(!weighted) effective_eb(0, ebs, ei);
    double max_value = 0, max_act = 0;
    size_t max_index = 0;
    for(size_t i = 0; i < n; i++){
        for(int k = 0; k < nv; k++) xi[vars[k]] = xptr[vars[k]][i];
        if(weighted) effective_eb(i, ebs, ei);
        bool masked = !mask_ptr || mask_ptr[i];
        double value;
        double est = Q::eval(xi, ei, masked, value);
        double act = fabs(value - qoi_ori[i]);
        if(act > max_act) max_act = act;
        if(max_value < est){
            max_value = est;
            max_index = i;
        }
    }
    max_act_error = max_act;
    max_est_error = max_value;
    if(!(max_value > tau)) return true;

    // tighten at the worst point
    const size_t i = max_index;
    for(int k = 0; k < nv; k++) xi[vars[k]] = xptr[vars[k]][i];
    const bool masked = !mask_ptr || mask_ptr[i];
    auto estimate = [&](const std::vector<double>& e){
        double value;
        effective_eb(i, e, ei);
        return Q::eval(xi, ei, masked, value);
    };
    std::vector<double> e = ebs;
    double est = max_value;
    if(inverse == Inverse::UNIFORM){
        while(est > tau){
            for(int v : q.vars) e[v] /= 1.5;
            est = estimate(e);
        }
    } else if(inverse == Inverse::COORDINATE_MIN){
        // tighten the variable(s) whose tightening lowers the estimate the most;
        // stop after a tenfold reduction of any bound (the next iteration continues)
        const double epsilon = 1e-6;
        while(est > tau){
            double min_error = std::numeric_limits<double>::infinity();
            std::vector<double> trial(NUM_VARS, 0);
            for(int v : q.vars){
                std::vector<double> t = e;
                t[v] /= 1.5;
                trial[v] = estimate(t);
                if(trial[v] < min_error) min_error = trial[v];
            }
            // an infinite estimate (a bound whose precondition fails) tightens every variable
            for(int v : q.vars){
                if(!std::isfinite(min_error) || fabs(trial[v] - min_error) < epsilon) e[v] /= 1.5;
            }
            est = min_error;
            bool stop = false;
            for(int v : q.vars) if(ebs[v] / e[v] > 10) stop = true;
            if(stop) break;
        }
    } else {
        // smooth proportional update: each variable is tightened in proportion to
        // how much its own tightening would lower the estimate (HPDC'26)
        const double factor_base = 1.5;
        const double alpha = 1.0 - 1.0 / factor_base;
        while(est > tau){
            std::vector<double> trial(NUM_VARS, 0);
            double sum_err = q.vars.size() * est;
            for(int v : q.vars){
                std::vector<double> t = e;
                t[v] /= 1.5;
                trial[v] = estimate(t);
                sum_err -= trial[v];
            }
            for(int v : q.vars){
                double w = (est - trial[v]) / sum_err;
                // an infinite estimate makes the proportional update undefined, and a zero
                // spread gives no direction: tighten every variable by the full factor instead
                if(!std::isfinite(sum_err) || !(sum_err > 0)) w = 1;
                e[v] *= (1.0 - alpha * w);
            }
            est = estimate(e);
        }
    }
    for(int v : q.vars) ebs[v] = e[v];
    return false;
}

template <class T>
bool hand_iteration(const HandQoI& q, const std::vector<const T*>& x, size_t n, const std::vector<unsigned char>& mask,
                    double tau, std::vector<double>& ebs, const std::vector<std::vector<int>>& weights, Inverse inverse,
                    const std::vector<double>& qoi_ori, double& max_act_error, double& max_est_error){
    return visit(q, [&](auto tag){
        return hand_iteration_impl<T>(tag, q, x, n, mask, tau, ebs, weights, inverse, qoi_ori, max_act_error, max_est_error);
    });
}

}

#endif
