#include "experiments.h"
#include "physics.h"
#include "mesh.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>

// TODO Refactor and clean these codes!
static Vec getEmptyLine(const LineModel& model, const SimConfig& cfg) {
   const double tol = Params{}.newton_tol;
   Vec q0 = model.static_solver_initial_guess();
   if (!model.static_solver(q0, tol))
      std::cerr << "warning: unloaded static solver did not converge\n";
   return q0;
}

// Max tension [N] of the unloaded line after pulling `pull` [m] of webbing at one anchor
static double unloadedMaxTension(const SimConfig& cfg, const Line& line, const double pull, const bool pull_right) {
   const Line pulled = pull_right ? pull_webbing(line, 0.0, pull) : pull_webbing(line, pull, 0.0);
   const DiscreteLine dline = discretize(pulled, elements_per_segment(pulled, cfg.L / (cfg.N - 1)));
   const LineModel model(dline, Vec2(0.0, 0.0), Vec2(cfg.L, 0.0), Params{});
   return model.max_tension_flat(getEmptyLine(model, cfg));
}

static Vec getInitialPosition(LineModel& model, const SimConfig& cfg, Eigen::Index& node_slackliner) {
   
   // Get initial state
   // TODO; Move newton tol
   const double tol = Params{}.newton_tol;
   Vec q0 = getEmptyLine(model, cfg);
   
   node_slackliner = model.nearest_node(q0, cfg.slackliner.x_coor);
   model.place_slackliner(node_slackliner, cfg.slackliner.m);
   if (!model.static_solver(q0, tol))
      std::cerr << "warning: loaded static solver did not converge\n";

   return q0;
}

bool Experiments::anyBreak(const std::vector<bool>& break_segments) {
   return std::find(break_segments.begin(), break_segments.end(), true) != break_segments.end();
}

static SimData simulateFall(const SimConfig& cfg, const DiscreteLine& start_line, const DiscreteLine& fall_line, const double t_end, const double dt_out, SolverType solvertype, const double atol, const double rtol) {

   Eigen::Index node_slackliner;
   const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
   const Params params{};

   LineModel start_model(start_line, anchor_left, anchor_right, params);

   Vec q0 = getInitialPosition(start_model, cfg, node_slackliner);

   LineModel fall_model(fall_line, anchor_left, anchor_right, params);

   FallSimulator fallsim(fall_model, cfg.slackliner, node_slackliner, solvertype);

   Vec v0(q0.size());
   v0.setZero();

   return fallsim.simulate(t_end, dt_out, as_const_mat2x(q0), as_const_mat2x(v0), atol, rtol);
 }

SimData Experiments::simulateLeashFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, SolverType solvertype, const double atol, const double rtol) {
   return simulateFall(cfg, dline, dline, t_end, dt_out, solvertype, atol, rtol);
 }

SimData Experiments::simulateBackupFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, const std::vector<bool>& break_segments, SolverType solvertype, const double atol, const double rtol) {
   assert(anyBreak(break_segments));
   return simulateFall(cfg, dline, with_broken_mainline(dline, break_segments), t_end, dt_out, solvertype, atol, rtol);
 }

Vec Experiments::energy_solve(const SimConfig& cfg, const DiscreteLine& dline, const std::vector<bool>& break_segments) {

   Eigen::Index node_slackliner;
   const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
   const Params params{};

   LineModel init_model(dline, anchor_left, anchor_right, params);

   Vec q = getInitialPosition(init_model, cfg, node_slackliner);

   const double g = 9.81; // TODO; put 9.82 somewhere sensible
                          //
   LineModel model(anyBreak(break_segments)? with_broken_mainline(dline, break_segments) : dline, anchor_left, anchor_right, params);
   model.place_slackliner(node_slackliner, cfg.slackliner.m);

   const double energy_goal = model.gravitational_energy_flat(q) + model.elastic_energy_flat(q) + cfg.slackliner.m * g * (cfg.slackliner.l_leg + cfg.slackliner.l_leash) ;

   const double tol = Params{}.newton_tol;
   q = model.static_solver_initial_guess();
   model.static_solver(q, tol);
   if (!model.static_solver_energy(q, tol, energy_goal))
      std::cerr << "warning: loaded static solver energy did not converge\n";

   return q;
}
  
Line Experiments::pull_to_tension(const SimConfig& cfg, const Line& line, const double tension, const bool pull_right) {
   const Segment& end = pull_right ? line.segments().back() : line.segments().front();

   // Tension increases with pull. pull_webbing can't pull a whole end segment
   const double max_pull = (1.0 - 1e-6) * end.L_main;

   if (unloadedMaxTension(cfg, line, 0.0, pull_right) >= tension)
      return line;

   // Bracket by doubling from a small pull, so the static solves stay at realistic pulls
   double lo = 0.0;
   double hi = std::min(0.1, max_pull);
   while (unloadedMaxTension(cfg, line, hi, pull_right) < tension) {
      if (hi == max_pull)
         throw std::invalid_argument("pull_to_tension: tension not reachable by pulling the end segment");
      lo = hi;
      hi = std::min(2.0 * hi, max_pull);
   }

   const double pull_tol = 1e-6; // [m]
   for (int it = 0; it < 100 && hi - lo > pull_tol; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (unloadedMaxTension(cfg, line, mid, pull_right) < tension) lo = mid;
      else hi = mid;
   }

   return pull_right ? pull_webbing(line, 0.0, hi) : pull_webbing(line, hi, 0.0);
}
