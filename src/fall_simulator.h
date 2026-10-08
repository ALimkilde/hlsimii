# pragma once

#include <vector>
#include "eigen_include.h"
#include "free_fall.h"
#include "physics.h"
#include "config.h"

struct SimData {
   double peak_force_leash;
   double peak_force_webbing;
   double peak_force_webbing_at_leash_fall; // /= peak_force_webbing => report shockload
   std::array<double, 2> peak_force_anchors;

   double lowest_point_slackliner;
   double lowest_point_webbing;

   // Some kind of logging of state of line + slackliner to plot later
};

class FallSimulator {
   public:
      explicit FallSimulator(const LineModel& linemodel, const FreeFall& freefall, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner);

      // Does time integration, detects events and restarts;
      // t will be overwritten with new time, q new positions and v new velocities
      SimData simulate(double& t, const double tend, Ref<Mat2X> q, Ref<Mat2X> v);

      // ====== Logging ====== 
      // Called once per time step
      // Returns False if we should stop..
      bool call_back(const double t, CRef<Mat2X> q, CRef<Mat2X> v, CRef<Mat2X> forces) const;

      // ====== Events ====== 

      // Detects when leash gets taut by checking freefall_ state y-only distance to node position
      // Event happens when sign of the function changes
      double event_leash_taut(const double t, const double y_slackliner) const;
      double event_leash_taut_flat(const double t, CRef<Vec> z) const;

      // Modifies velocities to preserve momentum
      void inelastic_collision(Ref<Mat2X> v);
      void inelastic_collision_flat(Ref<Vec> z);

      // Detects when to go back to free fall
      // Event happens when sign of the function changes
      double event_leash_loose(const double y_force_leash) const;
      double event_leash_loose_flat(const double t, CRef<Vec> z) const;

      // Free's the slackliner from the line and start free fall
      void free_slackliner(CRef<Vec2> q_slackliner, CRef<Vec2> v_slackliner);
      void free_slackliner_flat(CRef<Vec> z); 




   private:
      // TODO; what should be const or not (members and methods)
      LineModel linemodel_;
      FreeFall freefall_;
      const SlacklinerConfig slackliner_;
      Eigen::Index node_slackliner_;


};
