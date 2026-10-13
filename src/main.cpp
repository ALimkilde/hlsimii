#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

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
void write_fall_line_csv(const std::string& path, const SimData& data, Vec2 anchor_left, Vec2 anchor_right)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);

    out << "t,node,x,y\n";
    for (size_t f = 0; f < data.time.size(); ++f) {
        const Mat2X& q = data.position[f];
        const Eigen::Index num_nodes = q.cols();
        for (Eigen::Index i = -1; i <= num_nodes; ++i) {
            const Vec2 p = (i == -1) ? anchor_left : (i == num_nodes) ? anchor_right : Vec2(q.col(i));
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
        const Vec2 ring = data.position[f].col(data.node_slackliner);
        const double y_harness = data.y_slackliner[f];
        const bool taut = ring.y() - y_harness >= l_leash * (1.0 - 1e-9);
        out << data.time[f] << "," << ring.x() << "," << ring.y() << "," << y_harness << "," << taut << "\n";
    }
}

void print_summary(const std::string& label, const SimData& data)
{
    std::cout << label << " simulated to t = " << data.time.back() << " s, " << data.time.size() << " frames\n"
              << "  peak webbing tension " << data.peak_force_webbing / 1000.0 << " kN\n"
              << "  peak anchor forces   " << data.peak_force_anchors[0] / 1000.0 << " / "
              << data.peak_force_anchors[1] / 1000.0 << " kN\n"
              << "  lowest webbing       " << data.lowest_point_webbing << " m\n"
              << "  lowest slackliner    " << data.lowest_point_slackliner << " m\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: " << argv[0] << " <config.yaml> [webbings.yaml]\n";
        return 2;
    }
    const std::string webbings_path = argc == 3 ? argv[2] : "data/webbings.yaml";

    try {
        SimConfig cfg = parse_config(read_file(argv[1]));
        const WebbingCatalog webbings = parse_webbings(read_file(webbings_path));

        if (!cfg.pull_webbing) {
            std::cout << "Standing tension " << *cfg.tension_kN
                      << " kN (tension solve not implemented yet)\n";
            return 0;
        }

        const Line line = build_line(cfg, webbings);
        const double pull = *cfg.pull_webbing;
        const Line pulled = cfg.pull_side == "left" ? pull_webbing(line, pull, 0.0)
                                                    : pull_webbing(line, 0.0, pull);
        DiscreteLine dline = discretize(pulled, elements_per_segment(pulled, cfg.L / (cfg.N - 1)));

        const double dt_out = 2.0 / 60.0;   // frame rate of the written animation
        std::vector<bool> break_segments;
        for (const SegmentConfig& seg : cfg.segments) break_segments.push_back(seg.break_mainline);

        const FallResults results = Experiments::simulate(cfg, dline, cfg.T, dt_out, break_segments);

        // Must match the anchors used in Experiments::simulate
        const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
        const double l_leash = cfg.slackliner.l_leash;

        print_summary("Leash fall", results.leash);
        write_fall_line_csv("fall_line.csv", results.leash, anchor_left, anchor_right);
        write_fall_slackliner_csv("fall_slackliner.csv", results.leash, l_leash);
        std::cout << "Wrote fall_line.csv and fall_slackliner.csv\n";

        if (results.has_backup_fall) {
            print_summary("Backup fall", results.backup);
            write_fall_line_csv("fall_line_backup.csv", results.backup, anchor_left, anchor_right);
            write_fall_slackliner_csv("fall_slackliner_backup.csv", results.backup, l_leash);
            std::cout << "Wrote fall_line_backup.csv and fall_slackliner_backup.csv\n";
        } else {
            // Remove the files of an earlier run so the animation doesn't show a stale backup fall
            std::filesystem::remove("fall_line_backup.csv");
            std::filesystem::remove("fall_slackliner_backup.csv");
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
