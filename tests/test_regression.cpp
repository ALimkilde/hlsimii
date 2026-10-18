// Regression test: every case run with rk45 and rodas4 at the default tolerance must stay
// close to the golden written by make_golden, measured as a normalized RMS per column.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <chrono>
#include <cstdio>
#include <iostream>

#include "regression_case.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

const double test_tol = 1e-5;   // the solvers' default atol and rtol
const double rms_threshold = 1e-4;  // worst seen at test_tol: 1.2e-5 (example_backup, rk45, E_kin)

// Reads a header line of column names followed by rows of numbers. Lines starting with # are skipped.
regression::Trace read_csv(const std::string& path)
{
    std::istringstream file(regression::read_file(path));
    regression::Trace tr;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string field;
        if (tr.empty()) {
            while (std::getline(fields, field, ',')) tr.push_back({field, {}});
            continue;
        }
        for (regression::Column& col : tr) {
            if (!std::getline(fields, field, ',')) throw std::runtime_error(path + ": short row: " + line);
            col.values.push_back(std::stod(field));
        }
    }
    return tr;
}

struct Worst {
    std::string column;
    double rms = 0.0;
};

// Checks every column against the golden and returns the one furthest off
Worst check_against_golden(const std::string& fall, const regression::Trace& golden, const regression::Trace& run,
                           const char* solver)
{
    Worst worst;
    CHECK(run.size() == golden.size());
    CHECK(run[0].values.size() == golden[0].values.size());
    if (run.size() != golden.size() || run[0].values.size() != golden[0].values.size()) return worst;

    // Same output grid; a shifted grid would show up as a large RMS for a harmless reason
    CHECK(regression::normalized_rms(run[0].values, golden[0].values) < 1e-12);

    for (size_t c = 1; c < golden.size(); ++c) {
        CHECK(run[c].name == golden[c].name);
        const double e = regression::normalized_rms(run[c].values, golden[c].values);
        if (e > worst.rms) worst = {golden[c].name, e};
        if (e >= rms_threshold) {
            std::cerr << fall << " " << solver << " " << golden[c].name << ": normalized RMS " << e
                      << " >= " << rms_threshold << "\n";
            ++failures;
        }
    }
    return worst;
}

} // namespace

int main()
{
    std::printf("\n%-24s %-7s %8s %9s %9s  %s\n", "fall", "solver", "time [s]", "accepted", "rejected",
                "worst normalized RMS");
    for (const std::string& name : regression::cases) {
        for (const SolverType solver : {SolverType::RungeKutta45, SolverType::Rosenbrock43}) {
            const char* solver_name = solver == SolverType::RungeKutta45 ? "rk45" : "rodas4";

            const auto start = std::chrono::steady_clock::now();
            const FallResults r = regression::run_case(name, solver, test_tol);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

            // The time covers the whole case (both falls), so it is printed on the first row only
            bool first = true;
            for (const auto& [fall, trace] : regression::falls_of(name, r)) {
                const SimData& data = first ? r.leash : r.backup;
                const regression::Trace golden = read_csv(regression::data_dir + "golden/" + fall + ".csv");
                const Worst worst = check_against_golden(fall, golden, trace, solver_name);
                std::printf("%-24s %-7s %8s %9d %9d  %.1e (%s)\n", fall.c_str(), solver_name,
                            first ? std::to_string(seconds).substr(0, 5).c_str() : "",
                            data.n_accepts, data.n_rejects, worst.rms, worst.column.c_str());
                first = false;
            }
        }
    }
    std::printf("threshold %.0e at atol = rtol = %.0e\n\n", rms_threshold, test_tol);

    if (failures == 0) std::cout << "All tests passed\n";
    return failures == 0 ? 0 : 1;
}
