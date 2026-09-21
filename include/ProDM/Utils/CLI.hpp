#ifndef PRODM_UTILS_CLI_HPP
#define PRODM_UTILS_CLI_HPP

// Utils/CLI.hpp - shared pieces of the prodm_refactor / prodm_retrieve command line tools:
// option parsing, the method / approximator / encoder vocabularies, variable
// file naming, and the small helpers both tools need.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <climits>

#include "ProDM/Utils/IOUtils.hpp"
#include "ProDM/Utils/StatUtils.hpp"

#include "ProDM/Namespace.hpp"

namespace ProDM::CLI {

inline void fail(const std::string& msg){
    std::cerr << "error: " << msg << std::endl;
    exit(-1);
}

inline double now_seconds(){
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}

// ---------------------------------------------------------------------------
// option parsing: "--name" flags, "--name value" and "--name v1 v2 ..." options
// ---------------------------------------------------------------------------
class Options {
public:
    Options(int argc, char** argv, int num_positional){
        int i = 1;
        while(i < argc && positional.size() < (size_t)num_positional){
            if(std::string(argv[i]).rfind("--", 0) == 0) break;
            positional.push_back(argv[i++]);
        }
        std::string current;
        for(; i < argc; i++){
            std::string a = argv[i];
            if(a.rfind("--", 0) == 0){
                current = a.substr(2);
                if(values.count(current)) fail("option --" + current + " given twice");
                values[current] = {};
            } else if(current.empty()){
                fail("unexpected argument '" + a + "' (positional arguments come first)");
            } else {
                values[current].push_back(a);
            }
        }
    }
    bool has(const std::string& name) const { return values.count(name) > 0; }
    // consume-and-check: every recognized option calls one of the getters
    std::string get(const std::string& name, const std::string& def = "") {
        seen.insert(name);
        auto it = values.find(name);
        if(it == values.end()) return def;
        if(it->second.size() != 1) fail("option --" + name + " expects exactly one value");
        return it->second[0];
    }
    std::vector<std::string> get_list(const std::string& name) {
        seen.insert(name);
        auto it = values.find(name);
        if(it == values.end()) return {};
        std::vector<std::string> out;
        for(const auto& v : it->second){
            // allow both "--vars a,b,c" and "--vars a b c"
            size_t start = 0;
            while(true){
                size_t pos = v.find(',', start);
                out.push_back(v.substr(start, pos == std::string::npos ? std::string::npos : pos - start));
                if(pos == std::string::npos) break;
                start = pos + 1;
            }
        }
        return out;
    }
    bool flag(const std::string& name) {
        seen.insert(name);
        auto it = values.find(name);
        if(it == values.end()) return false;
        if(!it->second.empty()) fail("option --" + name + " takes no value");
        return true;
    }
    int get_int(const std::string& name, int def) {
        std::string v = get(name);
        return v.empty() ? def : atoi(v.c_str());
    }
    double get_double(const std::string& name, double def) {
        std::string v = get(name);
        return v.empty() ? def : atof(v.c_str());
    }
    std::string choice(const std::string& name, const std::vector<std::string>& allowed, const std::string& def) {
        std::string v = get(name, def);
        for(const auto& a : allowed) if(a == v) return v;
        std::string list;
        for(const auto& a : allowed) list += (list.empty() ? "" : "|") + a;
        fail("--" + name + " must be one of " + list + " (got '" + v + "')");
        return v;
    }
    void reject_unknown() const {
        for(const auto& kv : values){
            if(!seen.count(kv.first)) fail("unknown option --" + kv.first);
        }
    }
    std::vector<std::string> positional;
private:
    std::map<std::string, std::vector<std::string>> values;
    std::set<std::string> seen;
};

// ---------------------------------------------------------------------------
// vocabularies
// ---------------------------------------------------------------------------
enum class Method { MDR, PROAICD, PDR, PDR_DELTA };

inline Method parse_method(const std::string& s){
    if(s == "mdr") return Method::MDR;
    if(s == "proaicd") return Method::PROAICD;
    if(s == "pdr") return Method::PDR;
    if(s == "pdr-delta") return Method::PDR_DELTA;
    fail("--method must be one of mdr|proaicd|pdr|pdr-delta (got '" + s + "')");
    return Method::MDR;
}

// the approximators the command line tools compile in (the library also has Dummy, SZ2 and MGARD,
// reachable through the GE tools; each one adds its refactor and reconstructor instantiations)
enum class Approximator { SZ3, HPEZ, GE };

inline Approximator parse_approximator(const std::string& s){
    if(s == "sz3") return Approximator::SZ3;
    if(s == "hpez") return Approximator::HPEZ;
    if(s == "ge") return Approximator::GE;
    fail("--approximator must be one of sz3|hpez|ge (got '" + s + "')");
    return Approximator::HPEZ;
}

// ---------------------------------------------------------------------------
// variables: data lives in <data_dir>/<var><suffix>, the refactored streams of
// each variable in <refactor_dir>/<var>_refactored/
// ---------------------------------------------------------------------------
struct Variable {
    std::string name;
    std::string data_file;
    std::string refactor_dir;   // with trailing slash
};

inline std::vector<Variable> resolve_variables(const std::vector<std::string>& names, const std::string& data_dir,
                                               const std::string& refactor_dir, const std::string& suffix){
    if(names.empty()) fail("--vars is required");
    std::vector<Variable> vars;
    for(const auto& n : names){
        Variable v;
        v.name = n;
        // a variable may be given as a bare name (suffix appended) or as its file name
        struct stat st;
        std::string direct = data_dir + "/" + n;
        if(stat(direct.c_str(), &st) == 0 && !S_ISDIR(st.st_mode)){
            v.data_file = direct;
            size_t dot = n.find('.');
            v.name = (dot == std::string::npos) ? n : n.substr(0, dot);
        } else {
            v.data_file = data_dir + "/" + n + suffix;
        }
        v.refactor_dir = refactor_dir + "/" + v.name + "_refactored/";
        vars.push_back(v);
    }
    return vars;
}

inline std::string default_suffix(bool is_double){
    return is_double ? ".dat" : ".dat.f32";
}

template <class T>
std::vector<T> read_variable(const Variable& v, size_t expected_elements){
    size_t num = 0;
    auto data = ProDM::readfile<T>(v.data_file.c_str(), num);
    if(num == 0) fail("cannot read " + v.data_file);
    if(expected_elements && num != expected_elements){
        fail(v.data_file + " has " + std::to_string(num) + " elements, expected " + std::to_string(expected_elements));
    }
    return data;
}

inline std::vector<std::string> level_files(const std::string& dir, int num_levels){
    std::vector<std::string> files;
    for(int i = 0; i < num_levels; i++) files.push_back(dir + "level_" + std::to_string(i) + ".bin");
    return files;
}

// the number of levels stored by a refactor, read back from its metadata
// (multilevel layout: num_dims (1B), dims (4B each), num_levels (1B), ...;
//  approximation layout has an 8-byte approximator size between dims and num_levels)
inline int read_num_levels(const std::string& metadata_file, bool approximation_layout, int* num_dims_out = nullptr){
    size_t num_bytes = 0;
    auto metadata = ProDM::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
    if(num_bytes == 0) fail("cannot read " + metadata_file + "; run prodm_refactor first");
    int num_dims = metadata[0];
    size_t pos = 1 + num_dims * sizeof(uint32_t) + (approximation_layout ? sizeof(size_t) : 0);
    if(num_bytes <= pos) fail("metadata " + metadata_file + " is too short");
    if(num_dims_out) *num_dims_out = num_dims;
    return metadata[pos];
}

// proaicd metadata additionally stores the interpolation type after the level error bounds
template <class T>
int read_interpolation_type(const std::string& metadata_file){
    size_t num_bytes = 0;
    auto metadata = ProDM::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
    if(num_bytes == 0) fail("cannot read " + metadata_file + "; run prodm_refactor first");
    int num_dims = metadata[0];
    int num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
    size_t pos = num_dims * sizeof(uint32_t) + num_levels * sizeof(T) + 2;
    if(num_bytes <= pos) fail("metadata " + metadata_file + " is too short");
    return metadata[pos];
}

inline void make_dir(const std::string& path){
    mkdir(path.c_str(), 0777);
}

// absolute form of an existing path (the GE approximator locates the mesh block
// sizes from the refactor path, which must therefore be spelled out in full)
inline std::string absolute_path(const std::string& path){
    char buffer[PATH_MAX];
    if(realpath(path.c_str(), buffer) == nullptr) fail("cannot resolve path " + path);
    return std::string(buffer);
}

// the velocity mask of the SC'24 / HPDC'26 tools: nonzero total velocity
inline const std::vector<std::string>& velocity_names(){
    static const std::vector<std::string> names = {"VelocityX", "VelocityY", "VelocityZ"};
    return names;
}

inline bool is_velocity(const std::string& name){
    for(const auto& n : velocity_names()) if(n == name) return true;
    return false;
}

inline int find_variable(const std::vector<Variable>& vars, const std::string& name){
    for(size_t i = 0; i < vars.size(); i++) if(vars[i].name == name) return (int)i;
    return -1;
}

}

#endif
