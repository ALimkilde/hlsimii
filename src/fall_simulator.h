# pragma once

#include <array>
#include <limits>
#include <vector>
#include "eigen_include.h"
#include "free_fall.h"
#include "physics.h"
#include "config.h"

struct SimData {
   double peak_force_leash = 0.0;
   double peak_force_webbing = 0.0;
   double peak_force_webbing_at_leash_fall = 0.0; // /= peak_force_webbing => report shockload
   std::array<double, 2> peak_force_anchors = {0.0, 0.0};

   double lowest_point_slackliner = std::numeric_limits<double>::infinity();
   double lowest_point_webbing = std::numeric_limits<double>::infinity();

   // Frames at the output times: line state and slackliner (harness) height
   Eigen::Index node_slackliner = -1;
   std::vector<double> time;
   std::vector<Mat2X> position;
   std::vector<Mat2X> velocities;
   std::vector<Mat2X> net_forces;
   std::vector<double> y_slackliner;
   std::vector<double> gravitational_energy;
   std::vector<double> elastic_energy;
   std::vector<double> kinetic_energy;
};

class FallSimulator {
   public:
      explicit FallSimulator(const LineModel& linemodel, const SlacklinerConfig& slackliner, Eigen::Index node_slackliner);

      // Does time integration, detects events and restarts;
      // Frames are recorded every dt_out
      SimData simulate(const double tend, const double dt_out, CRef<Mat2X> q0, CRef<Mat2X> v0);

      // Harness height: free fall while the leash is slack, hanging l_leash below the node while taut
      double y_slackliner(const double t, CRef<Mat2X> q) const;

      // ====== Logging ====== 
      // Called once per output time.
      bool on_output(const double t, CRef<Vec> z, SimData& data, const bool attached) const;
      // Called once per time step
      bool on_step(const double t, CRef<Vec> z, SimData& data) const;

      // ====== Events ====== 

      // Detects when leash gets taut by checking freefall_ state y-only distance to node position
      // Event happens when sign of the function changes
      double event_leash_taut(const double t, CRef<Mat2X> q) const;
      double event_leash_taut_flat(const double t, CRef<Vec> z) const;

      // Modifies velocities to preserve momentum
      void inelastic_collision(double const t, Ref<Mat2X> v);

      // Detects when to go back to free fall
      // Event happens when sign of the function changes
      double event_leash_loose(const double y_acc_leash) const;
      double event_leash_loose_flat(const double t, CRef<Vec> z) const;

      // Free's the slackliner from the line and start free fall
      void free_slackliner(const double t, CRef<Mat2X> q, CRef<Mat2X> v);

   private:
      LineModel linemodel_;
      FreeFall freefall_;
      const SlacklinerConfig slackliner_;
      Eigen::Index node_slackliner_;


};
