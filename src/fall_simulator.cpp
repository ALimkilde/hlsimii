#include "fall_simulator.h"
#include "rk45.h"
#include "progress_bar.h"
#include <iostream>

FallSimulator::FallSimulator(const LineModel& linemodel, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner)
   : linemodel_(linemodel), 
     freefall_(FreeFall(Params{}.gravity.y(), 0.0, 0.0, 0.0)), 
     slackliner_(slackliner), 
     node_slackliner_(node_slackliner)
{
}

// TODO; warn if leash is lower than l_leg!
SimData FallSimulator::simulate(const double tend, const double dt_out, CRef<Mat2X> q0, CRef<Mat2X> v0)
{
  // Initialize
  SimData data;
  bool attached = false;

  // Slackliner leaves the line: the line starts at rest in the loaded shape,
  // the slackliner starts at rest with the harness l_leg above the standing node
  linemodel_.remove_slackliner();
  const double y0 = q0(1, node_slackliner_) + slackliner_.l_leg;
  const double t0 = 0.0;
  freefall_.set_initial_conditions(t0, y0, 0.0);

  // z = [q | v]
  Vec z(2 * q0.size());
  StateView state = state_view(z);
  state.q = q0;
  state.v = v0;

  data.node_slackliner = node_slackliner_;

  RKOptions opt;
  opt.dt_out = dt_out;

  ProgressBar bar;

  auto rhs    = [&](double t, CRef<Vec> z, Ref<Vec> dz) { linemodel_.rhs(t, z, dz); };
  auto output_obs = [&](double t, CRef<Vec> z) { bar.update((t - t0) / (tend - t0)); return on_output(t, z, data, attached); };
  auto step_obs   = [&](double t, CRef<Vec> z) { return on_step(t, z, data); };
  auto events = [&](double t, CRef<Vec> z, Ref<Vec> g) {
     g[0] = attached ? event_leash_loose_flat(t, z) : event_leash_taut_flat(t, z);
  };

  double t = t0;

  // Record the initial state here
  on_output(t, z, data, attached);

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

bool FallSimulator::on_output(const double t, CRef<Vec> z, SimData& data, const bool attached) const
{
   ConstStateView s = const_state_view(z);

   data.time.push_back(t);
   data.position.emplace_back(s.q);
   data.velocities.emplace_back(s.v);

   Mat2X F(2, s.q.cols());
   linemodel_.assemble(s.q, s.v, F);
   data.net_forces.emplace_back(F);

   data.y_slackliner.push_back(y_slackliner(t, s.q));
   data.elastic_energy.push_back( linemodel_.elastic_energy(s.q));

   double Ek = linemodel_.kinetic_energy(s.v);
   double Eg = linemodel_.gravitational_energy(s.q);
   if (!attached) {
      double ytmp = freefall_.get_y(t) + slackliner_.l_leash;
      Eg -= ytmp * freefall_.get_gy() * slackliner_.m;
      double vy = freefall_.get_vy(t);
      Ek += 0.5 * vy * vy * slackliner_.m;
   }
   data.gravitational_energy.push_back( Eg );
   data.kinetic_energy.push_back( Ek );

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
      const double T = linemodel_.edge_force(e, linemodel_.edge_vector(s.q, e), Vec2::Zero()).norm();
      data.peak_force_webbing = std::max(data.peak_force_webbing, T);
   }
   const double T_left  = linemodel_.edge_force(0, linemodel_.edge_vector(s.q, 0), Vec2::Zero()).norm();
   const double T_right = linemodel_.edge_force(num_nodes, linemodel_.edge_vector(s.q, num_nodes), Vec2::Zero()).norm();
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
   linemodel_.assemble(s.q, s.v, F); // Compute forces

   const double y_acc_leash = F(1, node_slackliner_) / linemodel_.get_mass(node_slackliner_);

   return event_leash_loose(y_acc_leash);

}

void FallSimulator::free_slackliner(const double t, CRef<Mat2X> q, CRef<Mat2X> v)
{
   linemodel_.remove_slackliner();

   double y_slackliner = q(1, node_slackliner_) - slackliner_.l_leash;
   const double nudge = 1e-12;

   freefall_.set_initial_conditions(t, y_slackliner+nudge, v(1, node_slackliner_));
}
