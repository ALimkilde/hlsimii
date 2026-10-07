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
    if (num_el < 2)
        throw std::invalid_argument("discretize: the line needs at least two elements (one free node)");

    DiscreteLine d;
    d.k_main.reserve(num_el);
    d.k_backup.reserve(num_el);
    d.l_main.reserve(num_el);
    d.l_backup.reserve(num_el);
    d.element_segment.reserve(num_el);

    // Lumped mass on all nodes 0..E, anchors included
    std::vector<double> lumped(num_el + 1, 0.0);

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
            lumped[e] += 0.5 * m;
            lumped[e + 1] += 0.5 * m;
        }
    }

    // Spread each anchor's mass evenly over the free nodes of its segment.
    // First segment: nodes 1..n_first. Last segment: nodes E-n_last..E-1.
    const int first_end = std::min(n_per_segment.front(), num_el - 1);
    for (int i = 1; i <= first_end; ++i)
        lumped[i] += lumped[0] / first_end;

    const int last_begin = std::max(num_el - n_per_segment.back(), 1);
    const int last_count = num_el - last_begin;
    for (int i = last_begin; i <= num_el - 1; ++i)
        lumped[i] += lumped[num_el] / last_count;

    // Keep only the free nodes 1..E-1
    d.node_mass.assign(lumped.begin() + 1, lumped.end() - 1);
    return d;
}
