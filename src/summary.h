#pragma once

#include <optional>
#include <ostream>

#include "experiments.h"

// Machine-readable version of the summary table main.cpp prints, for scripts.
// SI units throughout (m, N, kg); heights are y, so negative below the anchors.
//
// {
//   "setup_mass": kg,
//   "empty_line":   {"sag", "max_tension", "tension_left", "tension_right"},
//   "walking":      {"y_ring", "max_tension", "tension_left", "tension_right"},
//   "leash_fall":   {"y_harness_start", "y_harness_min",
//                    "peak_tension_left", "peak_tension_right", "peak_tension_leash"},
//   "backup_fall":  same keys as leash_fall plus
//                   "settled": {"y_harness", "tension_left", "tension_right"},
//                   or null when no mainline breaks
// }
void write_summary_json(std::ostream& out, const SimConfig& cfg, const ExperimentOutput& leash,
                        const std::optional<ExperimentOutput>& backup);
