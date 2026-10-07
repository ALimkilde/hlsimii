# pragma once

#include "eigen_include.h"

class FallingSlackliner {
   public:
      explicit FallingSlackliner(const double mass, const double t0, 
         const Vec2& p0, const Vec2& v0);

      void set_initial_conditions(const double t0, const Vec2& p0, const Vec2& v0) { t0_ = t0; p0_ = p0; v0_ = v0; }

      void get_state(const double t, Vec2& p, Vec2& v) const {

         // Do some maths

      }

      double get_mass() const { return mass_; }

   private:
      const double mass_;
      double t0_; // Time of p0 and v0
      Vec2 p0_;   // Position at t0
      Vec2 v0_;   // Velocity at t0

}
