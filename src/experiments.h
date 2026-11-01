#pragma once

#include "eigen_include.h"
#include "fall_simulator.h"

struct StaticData {
   StaticData(LineModel& model, Vec& q_in) : q(q_in) {
      max_tension = model.max_tension_flat(q);
      anchor_tension = model.anchor_tension_flat(q);
      sag = model.sag_flat(q);
   }

   Vec q;
   AnchorTension anchor_tension;
   double max_tension;
   double sag;
};

struct Anchors {
   explicit Anchors(const SimConfig& cfg): left(0.0, cfg.anchor_height_left), right(cfg.L, cfg.anchor_height_right) {}

   Vec2 left;
   Vec2 right;
};

struct ExperimentOutput {
   SimData dynamics;
   StaticData empty_line; // Static solution for empty line
   StaticData walking_line; // Static solution for walking line
   StaticData settled_line; // Static solution for settled line
   Anchors anchors;
   std::vector<int> element_segment;
   double setup_mass;
};

struct RunOptions {
   double dt_out = 1.0 / 60.0;   // frame rate of the written animation
   double atol = 1e-7;
   double rtol = 1e-3;
   SolverType solvertype = SolverType::Rosenbrock43;
};


namespace Experiments {
   // Pulls webbing at one anchor until the unloaded line's max tension is `tension` [N]
   Line pull_to_tension(const SimConfig& cfg, const Line& line, const double tension, const bool pull_right);

   // Setup a discreteline for simulation
   DiscreteLine setupTensionedLine(const SimConfig& cfg, const WebbingCatalog& webbings);

   ExperimentOutput simulateLeashFall(const SimConfig& cfg, const WebbingCatalog& webbings,  const RunOptions& opt);

   ExperimentOutput simulateBackupFall(const SimConfig& cfg, const WebbingCatalog& webbings, const RunOptions& opt);

   Vec energy_solve(const SimConfig& cfg, const DiscreteLine& dline, const std::vector<bool>& break_segments);

}
