#pragma once

#include "eigen_include.h"
#include "fall_simulator.h"

namespace Experiments {

   SimData simulate_leash_fall(SimConfig& cfg, DiscreteLine& dline, const double t_end, const double dt_out);

}
