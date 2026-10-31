#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"
#include "config_yaml.h"
#include "experiments.h"
#include "line.h"
#include "mesh.h"

namespace {

std::string read_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// One row per node per frame, anchors included (node -1 and num_nodes).
void write_fall_line_csv(const std::string& path, const ExperimentOutput& output)
{
    const SimData& data = output.dynamics;
    const Anchors& anchors = output.anchors;

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);

    out << "t,node,x,y\n";
    for (size_t f = 0; f < data.time.size(); ++f) {
        const Eigen::Map<const Mat2X> q = as_const_mat2x(data.frame(f));
        const Eigen::Index num_nodes = q.cols();
        for (Eigen::Index i = -1; i <= num_nodes; ++i) {
            const Vec2 p = (i == -1) ? anchors.left : (i == num_nodes) ? anchors.right : Vec2(q.col(i));
            out << data.time[f] << "," << i << "," << p.x() << "," << p.y() << "\n";
        }
    }
}

// One row per frame: leash ring (the node the slackliner falls onto) and harness.
// taut = 1 when the harness hangs a full leash length below the ring.
void write_fall_slackliner_csv(const std::string& path, const SimData& data, double l_leash)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);

    out << "t,x_ring,y_ring,y_harness,taut\n";
    for (size_t f = 0; f < data.time.size(); ++f) {
        const Vec2 ring = as_const_mat2x(data.frame(f)).col(data.node_slackliner);
        const double y_harness = data.y_slackliner[f];
        const bool taut = ring.y() - y_harness >= l_leash * (1.0 - 1e-9);
        out << data.time[f] << "," << ring.x() << "," << ring.y() << "," << y_harness << "," << taut << "\n";
    }
}

std::string two_decimals(double v)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << v;
    return ss.str();
}

std::string kN(double N) { return two_decimals(N / 1000.0); }

// Summary table: one row per line state, like the ISA backup fall calculator.
// Height is the leash ring for "Walking" and the harness for the falls; a fall's height is
// followed by the static energy estimate of its lowest point in parentheses.
const int label_w = 34, height_w = 16, force_w = 12;

void print_table_header()
{
    std::cout << std::left << std::setw(label_w) << "" << std::right << std::setw(height_w) << "Height [m]"
              << std::setw(force_w) << "Left [kN]" << std::setw(force_w) << "Right [kN]"
              << std::setw(force_w) << "Leash [kN]" << "\n";
}

void print_table_row(const std::string& label, const std::string& height, AnchorTension anchors, const std::string& leash = "-")
{
    std::cout << std::left << std::setw(label_w) << label << std::right << std::setw(height_w) << height
              << std::setw(force_w) << kN(anchors.left) << std::setw(force_w) << kN(anchors.right)
              << std::setw(force_w) << leash << "\n";
}

// Peaks of a simulated fall; q_energy is energy_solve's static estimate of its lowest point,
// where the harness hangs a leash length below the lowest node
void print_fall_row(const std::string& label, const SimData& d, const Vec& q_energy, double l_leash)
{
    const double y_energy = as_const_mat2x(q_energy).row(1).minCoeff() - l_leash;
    const std::string height = two_decimals(d.lowest_point_slackliner) + " (" + two_decimals(y_energy) + ")";
    print_table_row(label, height, {d.peak_force_anchors[0], d.peak_force_anchors[1]}, kN(d.peak_force_leash));
}

SolverType parse_solver(const std::string& name)
{
    if (name == "rk45") return SolverType::RungeKutta45;
    if (name == "rodas4") return SolverType::Rosenbrock43;
    throw std::runtime_error("unknown solver '" + name + "', expected rodas4 or rk45");
}

} // namespace

int main(int argc, char** argv)
{
    const std::string usage = std::string("usage: ") + argv[0] + " <config.yaml> [webbings.yaml] [--solver rodas4|rk45 (default rodas4)]\n";

    // Positional arguments in order, --solver anywhere
    std::vector<std::string> positional;
    std::string solver_name = "rodas4";
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--solver") {
            if (i + 1 == argc) { std::cerr << usage; return 2; }
            solver_name = argv[++i];
        } else {
            positional.push_back(arg);
        }
    }
    if (positional.empty() || positional.size() > 2) {
        std::cerr << usage;
        return 2;
    }
    const std::string config_path = positional[0];
    const std::string webbings_path = positional.size() == 2 ? positional[1] : "data/webbings.yaml";

    try {
        const SolverType solvertype = parse_solver(solver_name);
        SimConfig cfg = parse_config(read_file(config_path));
        const WebbingCatalog webbings = parse_webbings(read_file(webbings_path));

        const Line line = build_line(cfg, webbings);
        const bool pull_right = cfg.pull_side == "right";
        Line pulled = line;
        if (cfg.pull_webbing) {
            const double pull = *cfg.pull_webbing;
            pulled = pull_right ? pull_webbing(line, 0.0, pull) : pull_webbing(line, pull, 0.0);
        } else {
            pulled = Experiments::pull_to_tension(cfg, line, *cfg.tension_kN * 1000.0, pull_right);
        }
        DiscreteLine dline = discretize(pulled, elements_per_segment(pulled, cfg.L / (cfg.N - 1)));

        const double dt_out = 1.0 / 60.0;   // frame rate of the written animation
        const double atol = 1e-7, rtol = 1e-3;
        std::vector<bool> break_segments;
        for (const SegmentConfig& seg : cfg.segments) break_segments.push_back(seg.break_mainline);

        const bool has_backup_fall = Experiments::anyBreak(break_segments);

        const double l_leash = cfg.slackliner.l_leash;

        // The falls draw progress bars, so run everything before printing the table
        const std::vector<bool> no_breaks(break_segments.size(), false);
        const ExperimentOutput leash = Experiments::simulateLeashFall(cfg, dline, cfg.T, dt_out, solvertype, atol, rtol);
        const Vec q_energy_leash = Experiments::energy_solve(cfg, dline, no_breaks);
        write_fall_line_csv("fall_line.csv", leash);
        write_fall_slackliner_csv("fall_slackliner.csv", leash.dynamics, l_leash);

        std::optional<ExperimentOutput> backup;
        Vec q_energy_backup;
        if (has_backup_fall) {
            backup = Experiments::simulateBackupFall(cfg, dline, cfg.T, dt_out, break_segments, solvertype, atol, rtol);
            q_energy_backup = Experiments::energy_solve(cfg, dline, break_segments);
            write_fall_line_csv("fall_line_backup.csv", *backup);
            write_fall_slackliner_csv("fall_slackliner_backup.csv", backup->dynamics, l_leash);
        } else {
            // Remove the files of an earlier run so the animation doesn't show a stale backup fall
            std::filesystem::remove("fall_line_backup.csv");
            std::filesystem::remove("fall_slackliner_backup.csv");
        }

        const Eigen::Index ring = leash.dynamics.node_slackliner;
        const auto ring_y = [ring](const StaticData& s) { return as_const_mat2x(s.q)(1, ring); };

        print_table_header();
        print_table_row("Standing, no one on line (sag)", two_decimals(leash.empty_line.sag), leash.empty_line.anchor_tension);
        print_table_row("Walking (ring)", two_decimals(ring_y(leash.walking_line)), leash.walking_line.anchor_tension);
        print_fall_row("Leash fall (harness)", leash.dynamics, q_energy_leash, l_leash);
        if (backup) {
            print_fall_row("Backup fall - impact (harness)", backup->dynamics, q_energy_backup, l_leash);
            print_table_row("Backup fall - settled (harness)", two_decimals(ring_y(backup->settled_line) - l_leash),
                            backup->settled_line.anchor_tension);
        } else {
            std::cout << "Backup fall: no mainline break in this setup\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
