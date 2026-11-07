#include "fall_simulator.h"
#include "rk45.h"
#include "rodas4.h"
#include "progress_bar.h"
#include <iostream>

FallSimulator::FallSimulator(const LineModel& linemodel, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner, SolverType solvertype)
   : linemodel_(linemodel), 
     freefall_(FreeFall(Params{}.gravity.y(), 0.0, 0.0, 0.0)), 
     slackliner_(slackliner), 
     node_slackliner_(node_slackliner),
     solvertype_(solvertype)
{
}

// TODO; warn if leash is lower than l_leg!
SimData FallSimulator::simulate(const double tend, const double dt_out, CRef<Mat2X> q0, CRef<Mat2X> v0, const double atol, const double rtol)
{
  // Initialize
  SimData data{};
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
  data.num_nodes = q0.cols();

  ProgressBar bar;

  auto rhs    = [&](double t, CRef<Vec> z, Ref<Vec> dz) { linemodel_.rhs(t, z, dz); };
  auto rhsjac    = [&](double t, CRef<Vec> z, Ref<Vec> dz, Jac jacobian) { linemodel_.rhs(t, z, dz, &jacobian.Kq, &jacobian.Kv); };
  auto output_obs = [&](double t, CRef<Vec> z) { bar.update((t - t0) / (tend - t0)); return on_output(t, z, data); };
  auto step_obs   = [&](double t, CRef<Vec> z) { return on_step(t, z, data); };
  auto events = [&](double t, CRef<Vec> z, Ref<Vec> g) {
     g[0] = attached ? event_leash_loose_flat(t, z) : event_leash_taut_flat(t, z);
  };

  double t = t0;
  bool reason_was_event;

  // Record the initial state here
  on_output(t, z, data);

  while (t < tend) {
     if (solvertype_ == SolverType::RungeKutta45) {
        rk45::RKOptions opt;
        opt.dt_out = dt_out;
        opt.atol = atol;
        opt.rtol = rtol;

        rk45::RKResult r = rk45::integrate_rk45(rhs, events, 1, output_obs, step_obs, t, tend, z, opt);
        t = r.t;
        reason_was_event = r.reason == rk45::StopReason::Event;
        data.n_accepts += r.n_accepts;
        data.n_rejects += r.n_rejects;
     } else if (solvertype_ == SolverType::Rosenbrock43) {
        rodas4::Rodas4Options opt;
        opt.dt_out = dt_out;
        opt.atol = atol;
        opt.rtol = rtol;

        rodas4::Rodas4Result r = rodas4::integrate_rodas4<RodasLineModelLinSys>(rhs, rhsjac, events, 1, output_obs, step_obs, t, tend, z, opt);
        t = r.t;
        reason_was_event = r.reason == rodas4::StopReason::Event;
        data.n_accepts += r.n_accepts;
        data.n_rejects += r.n_rejects;
     }

     StateView s = state_view(z);

     if (!reason_was_event) break;
     if (!attached) inelastic_collision(t, s.v); else free_slackliner(t, s.q, s.v);
     attached = !attached;
  }

  return data;
}

bool FallSimulator::on_output(const double t, CRef<Vec> z, SimData& data) const
{

   on_step(t, z, data);

   ConstStateView s = const_state_view(z);

   data.time.push_back(t);

   // create flat vector from frames x num_nodes
   data.positions.insert(data.positions.end(), s.q.data(), s.q.data() + s.q.size());

   // create flat vector from frames x num_edges
   const Vec active = linemodel_.backupActivated(s.q); 
   data.backup_activated.insert(data.backup_activated.end(), active.data(), active.data() + active.size());

   AnchorTension anchors = linemodel_.anchor_tension(s.q);
   data.force_anchor_left.push_back(anchors.left);
   data.force_anchor_right.push_back(anchors.right);

   data.max_tension.push_back(linemodel_.max_tension(s.q));
   data.leash_tension.push_back(linemodel_.max_leash_tension(s.q));


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

   const double T_left  = linemodel_.edge_force(0, linemodel_.edge_vector(s.q, 0), Vec2::Zero()).norm();
   const double T_right = linemodel_.edge_force(num_nodes, linemodel_.edge_vector(s.q, num_nodes), Vec2::Zero()).norm();
   data.peak_force_anchors[0] = std::max(data.peak_force_anchors[0], T_left);
   data.peak_force_anchors[1] = std::max(data.peak_force_anchors[1], T_right);

   data.peak_force_webbing = std::max(data.peak_force_webbing, linemodel_.max_tension(s.q));
   data.peak_force_leash = std::max(data.peak_force_leash, linemodel_.max_leash_tension(s.q));

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

   const double m_node = linemodel_.get_mass(node_slackliner_);
   const double m_slackliner = slackliner_.m;

   v_node(1) = (m_node * v_node(1) + m_slackliner * freefall_.get_vy(t)) / ( m_node + m_slackliner);

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
   const Vec2 F = linemodel_.node_force(s.q, s.v, node_slackliner_); // Compute force
   const double y_acc_leash = F(1) / linemodel_.get_mass(node_slackliner_);

   return event_leash_loose(y_acc_leash);
}

void FallSimulator::free_slackliner(const double t, CRef<Mat2X> q, CRef<Mat2X> v)
{
   linemodel_.remove_slackliner();

   double y_slackliner = q(1, node_slackliner_) - slackliner_.l_leash;
   const double nudge = 1e-12;

   freefall_.set_initial_conditions(t, y_slackliner+nudge, v(1, node_slackliner_));
}
