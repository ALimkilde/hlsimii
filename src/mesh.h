#pragma once

#include <vector>

#include "line.h"

// The computational model: a chain of E springs (elements) between E+1 point masses (nodes).
// Nodes 0 and E are the anchors; they are not stored. Free node i (0..E-2) is mesh node i+1.
struct DiscreteLine {
    std::vector<double> k_main;    // spring constant of the main line [N/m], size E
    std::vector<double> k_backup;  // spring constant of the backup [N/m], size E
    std::vector<double> l_main;    // rest length of the main line [m], size E
    std::vector<double> l_backup;  // rest length of the backup [m], size E
    std::vector<double> node_mass; // lumped mass of the free nodes [kg], size E-1
    std::vector<int> element_segment; // segment id of each element, size E
    bool break_mainline = false;   // backup fall: drop the main line spring term everywhere

    int num_elements() const;
    int num_nodes() const;         // free nodes, E-1
};

// Number of elements N(s) per segment such that each element's main rest length is close to h [m]
std::vector<int> elements_per_segment(const Line& line, double h);

// Subdivides each segment s of the line into n_per_segment[s] equal elements
DiscreteLine discretize(const Line& line, const std::vector<int>& n_per_segment);
