// Tests for write_summary_json, the machine-readable summary scripts read.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "experiments.h"
#include "summary.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

// Small two-segment rig, segment 0 breaks. Short T keeps it cheap.
SimConfig rig()
{
    SimConfig cfg;
    cfg.L = 12.0;
    cfg.N = 13;
    cfg.T = 0.05;
    cfg.pull_webbing = 0.3;
    cfg.slackliner = {80.0, 1.0, 1.2, 6.0};
    cfg.segments = {{"main", "backup", 6.1, 6.4, true}, {"main", "backup", 6.1, 6.4, false}};
    return cfg;
}

WebbingCatalog webbings()
{
    return {{"main", Webbing{5.0, 10.0, 50.0}}, {"backup", Webbing{8.0, 10.0, 30.0}}};
}

// Value of the first "key": after `from` in the text; NaN if the key is missing
double value_after(const std::string& json, const std::string& key, size_t from = 0)
{
    const std::string needle = "\"" + key + "\": ";
    const size_t at = json.find(needle, from);
    if (at == std::string::npos) return std::nan("");
    return std::strtod(json.c_str() + at + needle.size(), nullptr);
}

void test_with_backup_fall()
{
    const SimConfig cfg = rig();
    const ExperimentOutput leash = Experiments::simulateLeashFall(cfg, webbings(), RunOptions{});
    const std::optional<ExperimentOutput> backup = Experiments::simulateBackupFall(cfg, webbings(), RunOptions{});

    std::ostringstream ss;
    write_summary_json(ss, cfg, leash, backup);
    const std::string json = ss.str();

    // Values round-trip exactly (max_digits10)
    CHECK(value_after(json, "setup_mass") == leash.setup_mass);
    CHECK(value_after(json, "sag") == leash.empty_line.sag);
    CHECK(value_after(json, "tension_left", json.find("\"walking\"")) == leash.walking_line.anchor_tension.left);

    const size_t at_leash = json.find("\"leash_fall\"");
    CHECK(value_after(json, "y_harness_start", at_leash) == leash.dynamics.y_slackliner.front());
    CHECK(value_after(json, "peak_tension_leash", at_leash) == leash.dynamics.peak_force_leash);

    const size_t at_backup = json.find("\"backup_fall\"");
    CHECK(value_after(json, "peak_tension_right", at_backup) == backup->dynamics.peak_force_anchors[1]);
    CHECK(value_after(json, "y_harness_min", at_backup) == backup->dynamics.lowest_point_slackliner);

    // Settled harness: a leash length below the ring of the settled line
    const double ring_settled = as_const_mat2x(backup->settled_line.q)(1, leash.dynamics.node_slackliner);
    CHECK(value_after(json, "y_harness", json.find("\"settled\"")) == ring_settled - cfg.slackliner.l_leash);

    CHECK(json.front() == '{' && json.find("null") == std::string::npos);
}

void test_without_backup_fall()
{
    const SimConfig cfg = rig();
    const ExperimentOutput leash = Experiments::simulateLeashFall(cfg, webbings(), RunOptions{});

    std::ostringstream ss;
    write_summary_json(ss, cfg, leash, std::nullopt);
    CHECK(ss.str().find("\"backup_fall\": null") != std::string::npos);
}

} // namespace

int main()
{
    test_with_backup_fall();
    test_without_backup_fall();

    if (failures == 0) std::cout << "All summary tests passed\n";
    return failures == 0 ? 0 : 1;
}
