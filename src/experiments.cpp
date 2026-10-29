#include "experiments.h"
#include "physics.h"
#include "mesh.h"
#include <iostream>
#include <algorithm>

// TODO Refactor and clean these codes!

FallResults Experiments::simulate(const SimConfig& cfg, const DiscreteLine& dline, const double t_end, const double dt_out, const std::vector<bool>& break_segments, SolverType solvertype, const double atol, const double rtol) {

   const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
   const Params params{};
   
   LineModel model(dline, anchor_left, anchor_right, params);
   
   // Get initial state
   // TODO; Move newton tol
   const double tol = Params{}.newton_tol;
   Vec q0 = model.static_solver_initial_guess();
   if (!model.static_solver(q0, tol))
      std::cerr << "warning: unloaded static solver did not converge\n";
   
   Eigen::Index node_slackliner = model.nearest_node(q0, cfg.slackliner.x_coor);
   model.place_slackliner(node_slackliner, cfg.slackliner.m);
   if (!model.static_solver(q0, tol))
      std::cerr << "warning: loaded static solver did not converge\n";
   
   FallSimulator fallsim(model, cfg.slackliner, node_slackliner, solvertype);

   Vec v0(q0.size());
   v0.setZero();

   FallResults fallresults{};

   fallresults.leash = fallsim.simulate(t_end, dt_out, as_const_mat2x(q0), as_const_mat2x(v0), atol, rtol);

   const bool any_break = std::find(break_segments.begin(), break_segments.end(), true) != break_segments.end();
   if ( any_break ) {

      LineModel broken_model(with_broken_mainline(dline, break_segments), anchor_left, anchor_right, params);

      FallSimulator fallsim_broken(broken_model, cfg.slackliner, node_slackliner, solvertype);

      fallresults.backup = fallsim_broken.simulate(t_end, dt_out, as_const_mat2x(q0), as_const_mat2x(v0), atol, rtol);
      fallresults.has_backup_fall = true;
   }

   return fallresults;
 }

Vec Experiments::energy_solve(const SimConfig& cfg, const DiscreteLine& dline, const std::vector<bool>& break_segments) {

   const Vec2 anchor_left(0.0, 0.0), anchor_right(cfg.L, 0.0);
   const Params params{};

   LineModel InitModel(dline, anchor_left, anchor_right, params);

   const double tol = Params{}.newton_tol;
   Vec q = InitModel.static_solver_initial_guess();
   if (!InitModel.static_solver(q, tol))
      std::cerr << "warning: unloaded static solver did not converge\n";

   Eigen::Index node_slackliner = InitModel.nearest_node(q, cfg.slackliner.x_coor);
   InitModel.place_slackliner(node_slackliner, cfg.slackliner.m);
   if (!InitModel.static_solver(q, tol))
      std::cerr << "warning: loaded static solver did not converge\n";

   const double g = 9.81; // TODO; put 9.82 somewhere sensible
                          //
   const bool any_break = std::find(break_segments.begin(), break_segments.end(), true) != break_segments.end();
   LineModel model(any_break? with_broken_mainline(dline, break_segments) : dline, anchor_left, anchor_right, params);
   model.place_slackliner(node_slackliner, cfg.slackliner.m);

   const double energy_goal = model.gravitational_energy_flat(q) + model.elastic_energy_flat(q) + cfg.slackliner.m * g * (cfg.slackliner.l_leg + cfg.slackliner.l_leash) ;

   q = model.static_solver_initial_guess();
   model.static_solver(q, tol);
   if (!model.static_solver_energy(q, tol, energy_goal))
      std::cerr << "warning: loaded static solver energy did not converge\n";

   return q;
}
