#include "mesh.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

// DiscreteLine

int DiscreteLine::num_elements() const { return static_cast<int>(k_main.size()); }
int DiscreteLine::num_nodes() const { return static_cast<int>(node_mass.size()); }

// Meshing

std::vector<int> elements_per_segment(const Line& line, double h)
{
    if (h <= 0.0)
        throw std::invalid_argument("elements_per_segment: target size must be positive");

    std::vector<int> n;
    n.reserve(line.segments().size());
    for (const Segment& s : line.segments())
        n.push_back(std::max(1, static_cast<int>(std::lround(s.L_main / h))));
    return n;
}

DiscreteLine discretize(const Line& line, const std::vector<int>& n_per_segment)
{
    const std::vector<Segment>& segments = line.segments();
    if (n_per_segment.size() != segments.size())
        throw std::invalid_argument("discretize: n_per_segment must have one entry per segment");
    for (int n : n_per_segment)
        if (n < 1)
            throw std::invalid_argument("discretize: each segment needs at least one element");

    const int num_el = std::accumulate(n_per_segment.begin(), n_per_segment.end(), 0);

    DiscreteLine d;
    d.k_main.reserve(num_el);
    d.k_backup.reserve(num_el);
    d.l_main.reserve(num_el);
    d.l_backup.reserve(num_el);
    d.element_segment.reserve(num_el);
    d.node_mass.assign(num_el + 1, 0.0);

    int e = 0;
    for (std::size_t s = 0; s < segments.size(); ++s) {
        const Segment& seg = segments[s];
        const int n = n_per_segment[s];
        const double l_main = seg.L_main / n;
        const double l_backup = seg.L_backup / n;
        const double m = seg.mass() / n;

        for (int j = 0; j < n; ++j, ++e) {
            d.k_main.push_back(seg.ea_main / l_main);
            d.k_backup.push_back(seg.ea_backup / l_backup);
            d.l_main.push_back(l_main);
            d.l_backup.push_back(l_backup);
            d.element_segment.push_back(static_cast<int>(s));

            // Lump half of the element mass on each of its end nodes
            d.node_mass[e] += 0.5 * m;
            d.node_mass[e + 1] += 0.5 * m;
        }
    }
    return d;
}
