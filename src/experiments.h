#pragma once

#include "eigen_include.h"
#include "fall_simulator.h"

struct FallResults {
   SimData leash;
   SimData backup;
   bool has_backup_fall = false;
};

namespace Experiments {
   FallResults simulate(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, const std::vector<bool>& break_segments, SolverType solvertype, const double atol, const double rtol);

}
