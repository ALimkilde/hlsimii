// Tests for Experiments::energy_solve, the static estimate of the lowest point of a fall.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.
//
// energy_solve only returns q, so each test rebuilds the same LineModel and checks the
// physics directly:
//   - energy balance: E(q) = E(start) + m g (l_leg + l_leash), the start being the slackliner
//     standing on the intact line and E measured on the line the fall happens on
//   - every node is in equilibrium, except the y-force on the slackliner's node (the leash pull)
//   - the slackliner's node ends up below where it stood (the lower of the two solutions)

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "experiments.h"
#include "mesh.h"
#include "physics.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

// 12 m span, six edges in two segments, main slightly shorter than the span so it is tensioned.
// Five free nodes, the middle one (node 2) at x = 6.
DiscreteLine six_edge_line()
{
    DiscreteLine d;
    d.k_main = std::vector<double>(6, 5000.0);
    d.k_backup = std::vector<double>(6, 3000.0);
    d.l_main = std::vector<double>(6, 1.99);
    d.l_backup = std::vector<double>(6, 2.1);
    d.main_width = std::vector<double>(6, 0.025);
    d.backup_width = std::vector<double>(6, 0.025);
    d.node_mass = Vec::Constant(5, 0.1);
    d.element_segment = {0, 0, 0, 1, 1, 1};
    return d;
}

SimConfig config()
{
    SimConfig cfg;
    cfg.L = 12.0;
    cfg.slackliner.m = 80.0;
    cfg.slackliner.l_leg = 1.0;
    cfg.slackliner.l_leash = 1.2;   // differs from l_leg so mixing them up shows
    cfg.slackliner.x_coor = 6.0;
    return cfg;
}

const Eigen::Index node_slackliner = 2;

double energy(const LineModel& m, const Vec& q) { return m.gravitational_energy(as_const_mat2x(q)) + m.elastic_energy(as_const_mat2x(q)); }

LineModel model_with_slackliner(const DiscreteLine& line, const SimConfig& cfg)
{
    const LineModel m(line, {0.0, 0.0}, {cfg.L, 0.0}, Params{});
    m.place_slackliner(node_slackliner, cfg.slackliner.m);
    return m;
}

// Same start as energy_solve and simulate: the slackliner standing on the intact line
Vec start_shape(const SimConfig& cfg)
{
    const LineModel m(six_edge_line(), {0.0, 0.0}, {cfg.L, 0.0}, Params{});
    const double tol = Params{}.newton_tol;
    Vec q = m.static_solver_initial_guess();
    CHECK(m.static_solver(q, tol));
    m.place_slackliner(node_slackliner, cfg.slackliner.m);
    CHECK(m.static_solver(q, tol));
    return q;
}

// The fall happens on fall_line (broken or intact); its energy is measured on that line
// from the intact start shape
void check_fall_bottom(const DiscreteLine& fall_line, const std::vector<bool>& break_segments)
{
    const SimConfig cfg = config();
    const LineModel m = model_with_slackliner(fall_line, cfg);
    const Vec q_start = start_shape(cfg);

    const Vec q = Experiments::energy_solve(cfg, six_edge_line(), break_segments);
    CHECK(q.size() == q_start.size());

    // Energy balance, relative to the energy released by the fall
    const double g = Params{}.gravity.norm();
    const double released = cfg.slackliner.m * g * (cfg.slackliner.l_leg + cfg.slackliner.l_leash);
    const double E_goal = energy(m, q_start) + released;
    CHECK(std::abs(energy(m, q) - E_goal) <= 1e-6 * released);

    // Net forces vanish everywhere except the y-force on the slackliner's node, which must pull up
    // (the line holds the slackliner up through the leash)
    Vec r(q.size());
    m.residual(q, r);
    Mat2X F = as_mat2x(r);
    const double leash_pull = F(1, node_slackliner);
    F(1, node_slackliner) = 0.0;
    CHECK(F.cwiseAbs().maxCoeff() <= 1e-4);
    CHECK(leash_pull > 0.0);

    // Lower solution: the node is below where the slackliner stood
    const double y_start = as_const_mat2x(q_start)(1, node_slackliner);
    const double y_fall = as_const_mat2x(q)(1, node_slackliner);
    CHECK(y_fall < y_start);
}

void test_main_line_fall()
{
    check_fall_bottom(six_edge_line(), {false, false});
}

// Breaking segment 0 leaves only the backup there. The fall starts on the intact line, so the
// main line's elastic energy in segment 0 is lost when it breaks.
void test_backup_fall()
{
    const std::vector<bool> break_segments = {true, false};
    check_fall_bottom(with_broken_mainline(six_edge_line(), break_segments), break_segments);
}

// The backup is softer and longer, so falling on it must end lower than falling on the main line
void test_backup_fall_is_lower()
{
    const SimConfig cfg = config();
    const Vec q_main = Experiments::energy_solve(cfg, six_edge_line(), {false, false});
    const Vec q_backup = Experiments::energy_solve(cfg, six_edge_line(), {true, false});
    CHECK(as_const_mat2x(q_backup)(1, node_slackliner) < as_const_mat2x(q_main)(1, node_slackliner));
}

// --- pull_to_tension ---

// 12 m span with two 6.1 m segments, so the unloaded line is slack until ~0.2 m is pulled
Line two_segment_line()
{
    const Webbing main{5.0, 10.0, 50.0};
    const Webbing backup{8.0, 10.0, 30.0};
    return Line({Segment(main, backup, 6.1, 6.4), Segment(main, backup, 6.1, 6.4)});
}

SimConfig pull_config()
{
    SimConfig cfg = config();
    cfg.N = 25;
    return cfg;
}

// Max tension of the unloaded line, computed the same way pull_to_tension does
double unloaded_tension(const SimConfig& cfg, const Line& line)
{
    const DiscreteLine dline = discretize(line, elements_per_segment(line, cfg.L / (cfg.N - 1)));
    const LineModel model(dline, Vec2(0.0, 0.0), Vec2(cfg.L, 0.0), Params{});
    Vec q = model.static_solver_initial_guess();
    model.static_solver(q, Params{}.newton_tol);
    return model.max_tension_flat(q);
}

// Pull a known amount, measure the tension, and check pull_to_tension recovers the pull
void test_pull_to_tension_round_trip()
{
    const SimConfig cfg = pull_config();
    const Line line = two_segment_line();
    const double pull = 0.3;

    for (const bool pull_right : {true, false}) {
        const Line known = pull_right ? pull_webbing(line, 0.0, pull) : pull_webbing(line, pull, 0.0);
        const double tension = unloaded_tension(cfg, known);
        CHECK(tension > 1000.0);

        const Line found = Experiments::pull_to_tension(cfg, line, tension, pull_right);
        const double L_end = pull_right ? found.segments().back().L_main : found.segments().front().L_main;
        const double L_other = pull_right ? found.segments().front().L_main : found.segments().back().L_main;
        CHECK(std::abs((6.1 - L_end) - pull) < 1e-4);
        CHECK(L_other == 6.1);
        CHECK(std::abs(unloaded_tension(cfg, found) - tension) < 1e-3 * tension);
    }
}

// A tension the end segment can't provide must throw
void test_pull_to_tension_unreachable()
{
    bool threw = false;
    try {
        Experiments::pull_to_tension(pull_config(), two_segment_line(), 1e12, true);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
}

// --- setup mass ---

// Small two-segment rig, built through the config like a real run. Short T keeps the
// simulation cheap; the setup mass doesn't depend on the dynamics.
SimConfig mass_config()
{
    SimConfig cfg = config();
    cfg.N = 13;
    cfg.T = 0.05;
    cfg.pull_webbing = 0.3;
    cfg.pull_side = "right";
    cfg.segments = {{"main", "backup", 6.1, 6.4, true}, {"main", "backup", 6.1, 6.4, false}};
    return cfg;
}

WebbingCatalog mass_webbings()
{
    return {{"main", Webbing{5.0, 10.0, 50.0}}, {"backup", Webbing{8.0, 10.0, 30.0}}};
}

double line_mass(const Line& line)
{
    double m = 0.0;
    for (const Segment& s : line.segments()) m += s.mass();
    return m;
}

// The setup mass is the webbing in the span after tensioning: the pulled webbing has left
// the span, so it's less than the line in the bag. Same rig for leash and backup falls.
void test_setup_mass()
{
    const SimConfig cfg = mass_config();
    const WebbingCatalog webbings = mass_webbings();
    const Line bag = build_line(cfg, webbings);
    const double expected = line_mass(pull_webbing(bag, 0.0, *cfg.pull_webbing));
    CHECK(expected < line_mass(bag));

    const ExperimentOutput leash = Experiments::simulateLeashFall(cfg, webbings, RunOptions{});
    CHECK(std::abs(leash.setup_mass - expected) <= 1e-12 * expected);

    const ExperimentOutput backup = Experiments::simulateBackupFall(cfg, webbings, RunOptions{});
    CHECK(std::abs(backup.setup_mass - expected) <= 1e-12 * expected);
}

} // namespace

int main()
{
    test_main_line_fall();
    test_backup_fall();
    test_backup_fall_is_lower();
    test_pull_to_tension_round_trip();
    test_pull_to_tension_unreachable();
    test_setup_mass();

    if (failures == 0) std::cout << "All experiments tests passed\n";
    return failures == 0 ? 0 : 1;
}
