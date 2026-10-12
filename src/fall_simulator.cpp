#include "fall_simulator.h"

FallSimulator::FallSimulator(const LineModel& linemodel, const FreeFall& freefall, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner)
   : linemodel_(linemodel), 
     freefall_(freefall), 
     slackliner_(slackliner), 
     node_slackliner_(node_slackliner)
{
}

// TODO; warn if leash is lower than l_leg!
SimData FallSimulator::simulate(double& t, const double tend, Ref<Mat2X> q, Ref<Mat2X> v)
{
}

bool FallSimulator::call_back(const double t, CRef<Mat2X> q, CRef<Mat2X> v) const
{
   // Implement this
   return false; 
}

bool FallSimulator::call_back_flat(const double t, CRef<Vec> z) const
{
   ConstStateView s = const_state_view(z);
   return call_back(t, s.q, s.v);
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
