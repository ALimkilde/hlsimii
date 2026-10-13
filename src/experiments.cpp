#include "experiments.h"
#include "physics.h"
#include <iostream>


SimData Experiments::simulate_leash_fall(SimConfig& cfg, DiscreteLine& dline, const double t_end, const double dt_out) {

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
   
   FallSimulator fallsim(model, cfg.slackliner, node_slackliner);

   Vec v0(q0.size());
   v0.setZero();

   return fallsim.simulate(t_end, dt_out, as_const_mat2x(q0), as_const_mat2x(v0));

 }
