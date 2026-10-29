#pragma once

#include "eigen_include.h"
#include "fall_simulator.h"

namespace Experiments {
   bool anyBreak(const std::vector<bool>& break_segments);

   SimData simulateLeashFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, SolverType solvertype, const double atol, const double rtol);

   SimData simulateBackupFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, const std::vector<bool>& break_segments, SolverType solvertype, const double atol, const double rtol);

   Vec energy_solve(const SimConfig& cfg, const DiscreteLine& dline, const std::vector<bool>& break_segments);

}
