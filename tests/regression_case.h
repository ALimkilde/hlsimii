#pragma once

// Shared by make_golden and test_regression: run one fixture config the way main.cpp does
// and flatten the result into named columns (one value per output frame).

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "config.h"
#include "config_yaml.h"
#include "experiments.h"
#include "line.h"
#include "mesh.h"

namespace regression {

const std::string data_dir = std::string(HLSIM_SOURCE_DIR) + "/tests/data/regression/";
const std::vector<std::string> cases = {"short_and_static", "example", "1km"};
const double dt_out = 1.0 / 60.0;  // same frame rate as main.cpp

struct Column {
    std::string name;
    std::vector<double> values;
};
using Trace = std::vector<Column>;

inline std::string read_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// One simulated fall, named for the golden file: <case>_leash.csv, <case>_backup.csv
struct Fall {
    std::string name;
    SimData data;
};

// Runs both falls as main.cpp does; tol is used as both atol and rtol
inline std::vector<Fall> run_case(const std::string& name, SolverType solver, double tol)
{
    const SimConfig cfg = parse_config(read_file(data_dir + name + ".yaml"));
    const WebbingCatalog webbings = parse_webbings(read_file(data_dir + "webbings.yaml"));

    RunOptions opt;
    opt.dt_out = dt_out;
    opt.solvertype = solver;
    opt.atol = opt.rtol = tol;

    std::vector<Fall> falls = {{name + "_leash", Experiments::simulateLeashFall(cfg, webbings, opt).dynamics}};
    if (has_backup_fall(cfg))
        falls.push_back({name + "_backup", Experiments::simulateBackupFall(cfg, webbings, opt).dynamics});
    return falls;
}

// The slackliner trace; the same quantities as fall_slackliner.csv
inline Trace trace_of(const SimData& data)
{
    Trace tr = {{"t", {}}, {"x_ring", {}}, {"y_ring", {}}, {"y_harness", {}}};
    for (size_t f = 0; f < data.time.size(); ++f) {
        const auto ring = as_const_mat2x(data.frame(f)).col(data.node_slackliner);
        const double row[] = {data.time[f], ring.x(), ring.y(), data.y_slackliner[f]};
        for (size_t c = 0; c < tr.size(); ++c) tr[c].values.push_back(row[c]);
    }
    return tr;
}

// sqrt(mean((a - ref)^2)) / sqrt(mean(ref^2)): unit free, so one threshold fits every column
inline double normalized_rms(const std::vector<double>& a, const std::vector<double>& ref)
{
    double diff2 = 0.0, ref2 = 0.0;
    for (size_t i = 0; i < ref.size(); ++i) {
        diff2 += (a[i] - ref[i]) * (a[i] - ref[i]);
        ref2 += ref[i] * ref[i];
    }
    return std::sqrt(diff2 / ref2);
}

} // namespace regression
