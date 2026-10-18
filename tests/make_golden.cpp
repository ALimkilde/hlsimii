// Writes the regression goldens to tests/data/regression/golden/.
// Runs every case with rk45 and rodas4 at a tight tolerance and only writes a golden when the
// two solvers agree, so the golden is close to the true solution of the model.
// Not a ctest: run it by hand after an intended change to the physics or the fixtures.

#include <cstdio>
#include <filesystem>
#include <iostream>

#include "regression_case.h"

namespace {

const double golden_tol = 1e-9;
const double agreement_threshold = 1e-6;  // normalized RMS between rk45 and rodas4, per column; worst seen 5e-8

void write_trace(const std::string& path, const regression::Trace& tr)
{
    std::FILE* out = std::fopen(path.c_str(), "w");
    if (!out) throw std::runtime_error("cannot write " + path);
    std::fprintf(out, "# rk45 and rodas4 at atol = rtol = %g, written by make_golden\n", golden_tol);
    for (size_t c = 0; c < tr.size(); ++c) std::fprintf(out, "%s%s", c ? "," : "", tr[c].name.c_str());
    std::fprintf(out, "\n");
    // Full double precision, so rounding doesn't hide the differences the test measures
    for (size_t i = 0; i < tr[0].values.size(); ++i) {
        for (size_t c = 0; c < tr.size(); ++c) std::fprintf(out, "%s%.17g", c ? "," : "", tr[c].values[i]);
        std::fprintf(out, "\n");
    }
    std::fclose(out);
}

// Prints the per-column normalized RMS and returns the largest
double compare(const regression::Trace& a, const regression::Trace& b)
{
    if (a[0].values.size() != b[0].values.size()) {
        std::cout << "    frame counts differ: " << a[0].values.size() << " vs " << b[0].values.size() << "\n";
        return INFINITY;
    }
    double worst = 0.0;
    for (size_t c = 1; c < a.size(); ++c) {  // skip t
        const double e = regression::normalized_rms(a[c].values, b[c].values);
        std::printf("    %-10s %.3e\n", a[c].name.c_str(), e);
        worst = std::max(worst, e);
    }
    return worst;
}

} // namespace

int main()
{
    const std::string out_dir = regression::data_dir + "golden/";
    std::filesystem::create_directories(out_dir);
    bool all_agree = true;

    for (const std::string& name : regression::cases) {
        const auto rk = regression::falls_of(name, regression::run_case(name, SolverType::RungeKutta45, golden_tol));
        const auto ro = regression::falls_of(name, regression::run_case(name, SolverType::Rosenbrock43, golden_tol));

        for (size_t k = 0; k < rk.size(); ++k) {
            std::cout << rk[k].first << ": normalized RMS rk45 vs rodas4\n";
            if (compare(rk[k].second, ro[k].second) > agreement_threshold) {
                std::cout << "    DISAGREE, golden not written\n";
                all_agree = false;
                continue;
            }
            write_trace(out_dir + rk[k].first + ".csv", rk[k].second);
        }
    }
    return all_agree ? 0 : 1;
}
