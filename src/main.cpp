#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "config.h"
#include "config_yaml.h"
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

double total_mass(const Line& line)
{
    double m = 0.0;
    for (const Segment& s : line.segments()) m += s.mass();
    return m;
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
