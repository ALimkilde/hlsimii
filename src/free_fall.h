# pragma once

struct FallState {
   double y;  // y - position
   double vy; // y - velocity
};

class FreeFall {
   public:
      explicit FreeFall(const double gy, const double t0, const double y0, const double vy0)
         : gy_(gy), t0_(t0), y0_(y0), vy0_(vy0) {}

      void set_initial_conditions(const double t0, const double y0, const double vy0) { t0_ = t0; y0_ = y0; vy0_ = vy0; }

      FallState get_state(const double t) const {
         const double delta_t = t - t0_;
         const double vy = gy_ * delta_t + vy0_;
         const double y = 0.5 * gy_ * delta_t*delta_t + vy0_ * delta_t + y0_;
         return FallState{y, vy};
      }

   private:
      const double gy_;
      double t0_;   // Time of y0 and vy0
      double y0_;   // Y - position at t0
      double vy0_;  // Y - velocity at t0

};
