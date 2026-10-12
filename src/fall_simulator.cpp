#include "fall_simulator.h"
#include "rk45.h"
#include <iostream>

FallSimulator::FallSimulator(const LineModel& linemodel, const FreeFall& freefall, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner)
   : linemodel_(linemodel), 
     freefall_(freefall), 
     slackliner_(slackliner), 
     node_slackliner_(node_slackliner)
{
}

// TODO; warn if leash is lower than l_leg!
SimData FallSimulator::simulate(double& t, const double tend, const double dt_out, Ref<Mat2X> q, Ref<Mat2X> v)
{
  // Initialize
  SimData data;
  bool attached = false;

  // Get initial state
  // Static shape with the slackliner standing on the node nearest x_coor (picked from the unloaded shape)
  // TODO; Move newton tol
  // TODO; Move this code to experiment.cpp at some point
  const double tol = Params{}.newton_tol;
  Vec q0 = linemodel_.static_solver_initial_guess();
  if (!linemodel_.static_solver(q0, tol))
     std::cerr << "warning: unloaded static solver did not converge\n";

  node_slackliner_ = linemodel_.nearest_node(q0, slackliner_.x_coor);
  linemodel_.place_slackliner(node_slackliner_, slackliner_.m);
  if (!linemodel_.static_solver(q0, tol))
     std::cerr << "warning: loaded static solver did not converge\n";

  // Slackliner leaves the line: the line starts at rest in the loaded shape,
  // the slackliner starts at rest with the harness l_leg above the standing node
  linemodel_.remove_slackliner();
  const double y0 = as_const_mat2x(q0)(1, node_slackliner_) + slackliner_.l_leg;
  freefall_.set_initial_conditions(t, y0, 0.0);

  // z = [q | v]
  Vec z(2 * q0.size());
  z.head(q0.size()) = q0;
  z.tail(q0.size()).setZero();

  data.node_slackliner = node_slackliner_;

  RKOptions opt;
  opt.dt_out = dt_out;

  auto rhs    = [&](double t, CRef<Vec> z, Ref<Vec> dz) { linemodel_.rhs(t, z, dz); };
  auto output_obs = [&](double t, CRef<Vec> z) { return on_output(t, z, data); };
  auto step_obs   = [&](double t, CRef<Vec> z) { return on_step(t, z, data); };
  auto events = [&](double t, CRef<Vec> z, Ref<Vec> g) {
     g[0] = attached ? event_leash_loose_flat(t, z) : event_leash_taut_flat(t, z);
  };

  // Record the initial state here
  on_output(t, z, data);

  while (t < tend) {
     RKResult r = rk45(rhs, events, 1, output_obs, step_obs, t, tend, z, opt);
     t = r.t;

     StateView s = state_view(z);

     if (r.reason != StopReason::Event) break;
     if (!attached) inelastic_collision(t, s.v); else free_slackliner(t, s.q, s.v);
     attached = !attached;
  }

  return data;
}

bool FallSimulator::on_output(const double t, CRef<Vec> z, SimData& data) const
{
   ConstStateView s = const_state_view(z);

   data.time.push_back(t);
   data.position.emplace_back(s.q);
   data.velocities.emplace_back(s.v);
   data.y_slackliner.push_back(y_slackliner(t, s.q));

   return true; // Never stop on output
}

bool FallSimulator::on_step(const double t, CRef<Vec> z, SimData& data) const
{
   ConstStateView s = const_state_view(z);
   const Eigen::Index num_nodes = s.q.cols();

   // Lowest points
   data.lowest_point_webbing = std::min(data.lowest_point_webbing, s.q.row(1).minCoeff());

   data.lowest_point_slackliner = std::min(data.lowest_point_slackliner, y_slackliner(t, s.q));

   // TODO; Move this to physics
   for (Eigen::Index e = 0; e <= num_nodes; ++e) {
      const double T = linemodel_.edge_force(e, linemodel_.edge_vector(s.q, e)).norm();
      data.peak_force_webbing = std::max(data.peak_force_webbing, T);
   }
   const double T_left  = linemodel_.edge_force(0, linemodel_.edge_vector(s.q, 0)).norm();
   const double T_right = linemodel_.edge_force(num_nodes, linemodel_.edge_vector(s.q, num_nodes)).norm();
   data.peak_force_anchors[0] = std::max(data.peak_force_anchors[0], T_left);
   data.peak_force_anchors[1] = std::max(data.peak_force_anchors[1], T_right);

   // TODO peak_force_leash and peak_force_webbing_at_leash_fall: need to know if the slackliner is attached

   return true; // Never stop on step
}

double FallSimulator::y_slackliner(const double t, CRef<Mat2X> q) const
{
   // While free, the slackliner is above the taut-leash position; while attached, the stale
   // free fall parabola has fallen below it. So the max is the slackliner position in both phases.
   return std::max(freefall_.get_y(t), q(1, node_slackliner_) - slackliner_.l_leash);
}

double FallSimulator::event_leash_taut(const double t, CRef<Mat2X> q) const
{
   const double y_leash_ring = q(1, node_slackliner_);
   const double y_slackliner = freefall_.get_y(t);

   return y_leash_ring - slackliner_.l_leash - y_slackliner;
}

double FallSimulator::event_leash_taut_flat(const double t, CRef<Vec> z) const
{
   ConstStateView s = const_state_view(z);
   return event_leash_taut(t, s.q);
}

void FallSimulator::inelastic_collision(double const t, Ref<Mat2X> v)
{
   auto v_node = v.col(node_slackliner_); // Returns view

   Vec2 v_slackliner(0.0, freefall_.get_vy(t));

   const double m_node = linemodel_.get_mass(node_slackliner_);
   const double m_slackliner = slackliner_.m;

   v_node = (m_node * v_node + m_slackliner * v_slackliner) / ( m_node + m_slackliner);

   linemodel_.place_slackliner(node_slackliner_, m_slackliner);
}

double FallSimulator::event_leash_loose(const double y_acc_leash) const
{
   // Event happens when y_acc_leash < gy
   // Remember gy = -9.81
   return freefall_.get_gy() - y_acc_leash; 
}

double FallSimulator::event_leash_loose_flat(const double t, CRef<Vec> z) const
{
   ConstStateView s = const_state_view(z);
   const Eigen::Index N = s.q.cols();

   // TODO Can be optimized to only compute acceleration of specific node.
   Mat2X F(2, N);
   // TODO Add velocity here once drag and KV is added
   linemodel_.assemble(s.q, F); // Compute forces

   const double y_acc_leash = F(1, node_slackliner_) / linemodel_.get_mass(node_slackliner_);

   return event_leash_loose(y_acc_leash);

}

void FallSimulator::free_slackliner(const double t, CRef<Mat2X> q, CRef<Mat2X> v)
{
   linemodel_.remove_slackliner();

   double y_slackliner = q(1, node_slackliner_) - slackliner_.l_leash;

   freefall_.set_initial_conditions(t, y_slackliner , v(1, node_slackliner_));
}
