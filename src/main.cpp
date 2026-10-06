#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "config.h"
#include "config_yaml.h"
#include "line.h"
#include "mesh.h"
#include "physics.h"

namespace {

std::string read_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

double total_mass(const Line& line)
{
    double m = 0.0;
    for (const Segment& s : line.segments()) m += s.mass();
    return m;
}

// One row per node, anchors included. T [N] is the tension in the edge to the
// right of the node, so the right anchor gets nan. m_point [kg] is the point mass
// on the node: point_mass on free node point_node, 0 elsewhere (point_node = -1: none).
void write_shape_csv(const std::string& path, const LineModel& model, const Vec& q,
                     Vec2 anchor_left, Vec2 anchor_right,
                     int point_node = -1, double point_mass = 0.0)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);

    Eigen::Map<const Mat2X> Q(q.data(), 2, q.size() / 2);
    const Eigen::Index num_nodes = Q.cols();

    out << "x,y,T,m_point\n";
    for (Eigen::Index i = -1; i <= num_nodes; ++i) {
        const Vec2 p = (i == -1) ? anchor_left : (i == num_nodes) ? anchor_right : Vec2(Q.col(i));
        out << p.x() << "," << p.y() << ",";
        if (i < num_nodes) out << model.edge_force(i + 1, model.edge_vector(Q, i + 1)).norm();
        else out << "nan";
        out << "," << (i == point_node ? point_mass : 0.0) << "\n";
    }
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
        const SimConfig cfg = parse_config(read_file(argv[1]));
        const WebbingCatalog webbings = parse_webbings(read_file(webbings_path));
        const Line line = build_line(cfg, webbings);
        std::cout << "Line with " << line.num_segments() << " segments, total mass "
                  << total_mass(line) << " kg\n";

        if (cfg.pull_webbing) {
            const double pull = *cfg.pull_webbing;
            const Line pulled = cfg.pull_side == "left" ? pull_webbing(line, pull, 0.0)
                                                        : pull_webbing(line, 0.0, pull);
            const DiscreteLine d = discretize(pulled, elements_per_segment(pulled, cfg.L / (cfg.N - 1)));
            std::cout << "Pulled " << pull << " m on the " << cfg.pull_side << ": mesh with "
                      << d.num_elements() << " elements\n";

            const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
            const Params params{};
            LineModel model(d, anchor_left, anchor_right, params);
            Vec q = model.static_solver_initial_guess();
            if (!model.static_solver(q, params.newton_tol))
                std::cerr << "warning: static solver did not converge, writing last iterate\n";
            write_shape_csv("static_shape.csv", model, q, anchor_left, anchor_right);
            std::cout << "Wrote static_shape.csv\n";

            // Hang the slackliner on the node nearest x_coor, starting from the unloaded shape
            const int slackliner_node = model.nearest_node(q, cfg.slackliner.x_coor);
            std::cout << "Slackliner (" << cfg.slackliner.m << " kg) at node " << slackliner_node << "\n";
            model.place_slackliner(slackliner_node, cfg.slackliner.m);
            if (!model.static_solver(q, params.newton_tol))
                std::cerr << "warning: loaded static solver did not converge, writing last iterate\n";
            write_shape_csv("static_shape_loaded.csv", model, q, anchor_left, anchor_right,
                            slackliner_node, cfg.slackliner.m);
            std::cout << "Wrote static_shape_loaded.csv\n";
        } else {
            std::cout << "Standing tension " << *cfg.tension_kN
                      << " kN (tension solve not implemented yet)\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
