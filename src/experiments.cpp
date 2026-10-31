#include "experiments.h"
#include "physics.h"
#include "mesh.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>

static Params getParams(const SimConfig& cfg) 
{
   Params p;
   
   p.tau_KV *= cfg.kv_scale;
   p.c_drag *= cfg.drag_scale;

   return p;
}

static Vec getStaticPosition(Vec q0, const LineModel& model) {
   
   const double tol = Params{}.newton_tol;
   if (!model.static_solver(q0, tol))
      std::cerr << "warning: loaded static solver did not converge\n";

   return q0;
}

static Vec getEmptyLine(const LineModel& model) {
   Vec q0 = model.static_solver_initial_guess();
   return getStaticPosition(q0, model);
}

// Max tension [N] of the unloaded line after pulling `pull` [m] of webbing at one anchor
static double unloadedMaxTension(const SimConfig& cfg, const Line& line, const double pull, const bool pull_right) {
   const Line pulled = pull_right ? pull_webbing(line, 0.0, pull) : pull_webbing(line, pull, 0.0);
   const DiscreteLine dline = discretize(pulled, elements_per_segment(pulled, cfg.L / (cfg.N - 1)));
   const Anchors anchors(cfg);
   const LineModel model(dline, anchors.left, anchors.right, Params{});
   return model.max_tension_flat(getEmptyLine(model));
}

Eigen::Index placeSlackliner(Vec q0, LineModel& model, const SimConfig& cfg) 
{
   Eigen::Index node_slackliner = model.nearest_node(q0, cfg.slackliner.x_coor);
   model.place_slackliner(node_slackliner, cfg.slackliner.m);
   
   return node_slackliner;

}

bool Experiments::anyBreak(const std::vector<bool>& break_segments) {
   return std::find(break_segments.begin(), break_segments.end(), true) != break_segments.end();
}

static ExperimentOutput simulateFall(const SimConfig& cfg, const DiscreteLine& start_line, const DiscreteLine& fall_line, const double t_end, const double dt_out, SolverType solvertype, const double atol, const double rtol) {

   const Anchors anchors(cfg);
   const Params params = getParams(cfg);

   LineModel start_model(start_line, anchors.left, anchors.right, params);

   Vec q_empty = getEmptyLine(start_model);
   StaticData empty_data(start_model, q_empty);

   Eigen::Index node_slackliner = placeSlackliner(q_empty, start_model, cfg);

   Vec q_walking = getStaticPosition(q_empty, start_model);
   StaticData walking_data(start_model, q_walking);

   LineModel fall_model(fall_line, anchors.left, anchors.right, params);

   FallSimulator fallsim(fall_model, cfg.slackliner, node_slackliner, solvertype);

   Vec v0(q_walking.size());
   v0.setZero();

   SimData dynamics = fallsim.simulate(t_end, dt_out, as_const_mat2x(q_walking), as_const_mat2x(v0), atol, rtol);

   fall_model.place_slackliner(node_slackliner, cfg.slackliner.m);
   Vec q_settled = getStaticPosition(dynamics.lastFrame(), fall_model);
   StaticData settled_data(fall_model, q_settled);

   return {std::move(dynamics), empty_data, walking_data, settled_data, anchors};
 }

ExperimentOutput Experiments::simulateLeashFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, SolverType solvertype, const double atol, const double rtol) {
   return simulateFall(cfg, dline, dline, t_end, dt_out, solvertype, atol, rtol);
 }

ExperimentOutput Experiments::simulateBackupFall(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, const std::vector<bool>& break_segments, SolverType solvertype, const double atol, const double rtol) {
   assert(anyBreak(break_segments));
   return simulateFall(cfg, dline, with_broken_mainline(dline, break_segments), t_end, dt_out, solvertype, atol, rtol);
 }

Vec Experiments::energy_solve(const SimConfig& cfg, const DiscreteLine& dline, const std::vector<bool>& break_segments) {

   const Anchors anchors(cfg);
   const Params params{};

   LineModel init_model(dline, anchors.left, anchors.right, params);

   Vec q = getEmptyLine(init_model);
   Eigen::Index node_slackliner = placeSlackliner(q, init_model, cfg);
   q = getStaticPosition(q, init_model);

   const double g = 9.81; // TODO; put 9.82 somewhere sensible
                          //
   LineModel model(anyBreak(break_segments)? with_broken_mainline(dline, break_segments) : dline, anchors.left, anchors.right, params);
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
