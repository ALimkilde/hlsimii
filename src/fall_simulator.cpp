#include "fall_simulator.h"

FallSimulator::FallSimulator(const LineModel& linemodel, const FreeFall& freefall, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner)
   : linemodel_(linemodel), 
     freefall_(freefall), 
     slackliner_(slackliner), 
     node_slackliner_(node_slackliner)
{
}


SimData FallSimulator::simulate(double& t, const double tend, Ref<Mat2X> q, Ref<Mat2X> v)
{
}

bool FallSimulator::call_back(const double t, CRef<Mat2X> q, CRef<Mat2X> v, CRef<Mat2X> forces) const
{
}

double FallSimulator::event_leash_taut(const double t, const double y_slackliner) const
{
}

double FallSimulator::event_leash_taut_flat(const double t, CRef<Vec> z) const
{
}

void FallSimulator::inelastic_collision(Ref<Mat2X> v)
{
}

void FallSimulator::inelastic_collision_flat(Ref<Vec> z)
{
}

double FallSimulator::event_leash_loose(const double y_force_leash) const
{
}

double FallSimulator::event_leash_loose_flat(const double t, CRef<Vec> z) const
{
}

void FallSimulator::free_slackliner(CRef<Vec2> q_slackliner, CRef<Vec2> v_slackliner)
{
}

void FallSimulator::free_slackliner_flat(CRef<Vec> z)
{
}
