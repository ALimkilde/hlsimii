#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include "eigen_include.h"

// Adaptive explicit Runge-Kutta integration, structured like DifferentialEquations.jl:
//   - Stepper (DP5): one trial step + error estimate, accept, dense output. Method specific.
//   - PIController:  step size selection from the error estimate.
//   - integrate():   the generic loop. Step control, events on the interpolant, output grid
//                    and observers. Reusable by other steppers (e.g. RODAS4).
// rk45() is a convenience wrapper that runs integrate() with a DP5 stepper.
//

namespace rk45 {

struct RKOptions {
   double rtol = 1e-5, atol = 1e-5;
   double h0 = 1e-4, hmax = std::numeric_limits<double>::infinity();
   int max_steps = 1'000'000;
   double dt_out = 0.0;   // Output on the global grid t_k = k * dt_out; <= 0 disables output
};

enum class StopReason { ReachedEnd, Event, ObserverStop, StepTooSmall, MaxSteps };

struct RKResult {
   StopReason reason;
   double t;          // where integration stopped
   int event_index;   // which event fired (-1 if none)
};


// ====== Step size control ======

// PI controller (Julia's PIController). Error estimates are normalized: a step is accepted if err <= 1.
struct PIController {
   double beta1, beta2;            // Exponents on the current and previous error
   double qmin = 0.2, qmax = 10.0; // Bounds on the step size ratio h_new / h
   double gamma = 0.9;             // Safety factor
   double err_prev = 1e-4;         // Error of the last accepted step

   // Step size for the next step after accepting a step of size h with error err
   double accept(const double err, const double h) {
      const double q = (err == 0.0) ? 1.0 / qmax
                                    : std::clamp(std::pow(err, beta1) / std::pow(err_prev, beta2) / gamma,
                                                 1.0 / qmax, 1.0 / qmin);
      err_prev = std::max(err, 1e-4);
      return h / q;
   }

   // Step size to retry with after rejecting a step of size h with error err (may be NaN)
   double reject(const double err, const double h) const {
      if (!std::isfinite(err)) return h * qmin;
      return h / std::min(1.0 / qmin, std::pow(err, beta1) / gamma);
   }
};


// ====== Dormand-Prince 5(4) stepper ======

// 5th order solution, 4th order embedded error estimate, FSAL (the last stage is f at the new state)
// and a 4th order continuous extension. Same method as Julia's DP5 and Hairer's DOPRI5.
class DP5 {
   public:
      explicit DP5(const Eigen::Index n)
         : y_(n), y_new_(n), y_prev_(n), tmp_(n), err_(n),
           k1_(n), k2_(n), k3_(n), k4_(n), k5_(n), k6_(n), k7_(n),
           r2_(n), r3_(n), r4_(n), r5_(n) {}

      // Julia's DP5 controller defaults
      static PIController default_controller() { return PIController{1.0/5.0 - 0.75*0.04, 0.04}; }

      // Start at (t0, z0). Must be called after every discontinuity, since it resets the FSAL stage.
      template <class Rhs>
      void init(Rhs& rhs, const double t0, CRef<Vec> z0) {
         t_ = t_prev_ = t0;
         y_ = z0;
         y_prev_ = z0;
         rhs(t_, y_, k1_);
      }

      // Trial step of size h from the current state. Returns the normalized RMS error estimate.
      template <class Rhs>
      double try_step(Rhs& rhs, const double h, const double rtol, const double atol) {
         h_ = h;

         tmp_ = y_ + h * (a21*k1_);
         rhs(t_ + c2*h, tmp_, k2_);
         tmp_ = y_ + h * (a31*k1_ + a32*k2_);
         rhs(t_ + c3*h, tmp_, k3_);
         tmp_ = y_ + h * (a41*k1_ + a42*k2_ + a43*k3_);
         rhs(t_ + c4*h, tmp_, k4_);
         tmp_ = y_ + h * (a51*k1_ + a52*k2_ + a53*k3_ + a54*k4_);
         rhs(t_ + c5*h, tmp_, k5_);
         tmp_ = y_ + h * (a61*k1_ + a62*k2_ + a63*k3_ + a64*k4_ + a65*k5_);
         rhs(t_ + h, tmp_, k6_);
         y_new_ = y_ + h * (b1*k1_ + b3*k3_ + b4*k4_ + b5*k5_ + b6*k6_);
         rhs(t_ + h, y_new_, k7_);

         // Difference between the 5th and 4th order solutions, scaled per component
         err_ = h * (e1*k1_ + e3*k3_ + e4*k4_ + e5*k5_ + e6*k6_ + e7*k7_);
         err_.array() /= atol + rtol * y_.array().abs().max(y_new_.array().abs());
         return std::sqrt(err_.squaredNorm() / static_cast<double>(err_.size()));
      }

      // Accept the last trial step: advance t and the state, and build the dense output for it
      void accept() {
         r2_ = y_new_ - y_;
         r3_ = h_ * k1_ - r2_;
         r4_ = r2_ - h_ * k7_ - r3_;
         r5_ = h_ * (d1*k1_ + d3*k3_ + d4*k4_ + d5*k5_ + d6*k6_ + d7*k7_);

         t_prev_ = t_;
         t_ += h_;
         y_prev_.swap(y_);
         y_.swap(y_new_);
         k1_.swap(k7_);   // FSAL
      }

      // State at t in [t_prev(), t()], from the last accepted step
      void interpolate(const double t, Ref<Vec> out) const {
         const double th  = (t - t_prev_) / (t_ - t_prev_);
         const double th1 = 1.0 - th;
         out = y_prev_ + th * (r2_ + th1 * (r3_ + th * (r4_ + th1 * r5_)));
      }

      double t() const { return t_; }
      double t_prev() const { return t_prev_; }
      const Vec& state() const { return y_; }

   private:
      // Butcher tableau
      static constexpr double c2 = 1.0/5.0, c3 = 3.0/10.0, c4 = 4.0/5.0, c5 = 8.0/9.0;

      static constexpr double a21 = 1.0/5.0;
      static constexpr double a31 = 3.0/40.0,       a32 = 9.0/40.0;
      static constexpr double a41 = 44.0/45.0,      a42 = -56.0/15.0,      a43 = 32.0/9.0;
      static constexpr double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0,
                              a54 = -212.0/729.0;
      static constexpr double a61 = 9017.0/3168.0,  a62 = -355.0/33.0,     a63 = 46732.0/5247.0,
                              a64 = 49.0/176.0,     a65 = -5103.0/18656.0;

      // 5th order weights (b2 = 0)
      static constexpr double b1 = 35.0/384.0, b3 = 500.0/1113.0, b4 = 125.0/192.0,
                              b5 = -2187.0/6784.0, b6 = 11.0/84.0;

      // Error weights: 5th minus 4th order (e2 = 0)
      static constexpr double e1 = 71.0/57600.0, e3 = -71.0/16695.0, e4 = 71.0/1920.0,
                              e5 = -17253.0/339200.0, e6 = 22.0/525.0, e7 = -1.0/40.0;

      // Dense output (Hairer's DOPRI5, d2 = 0)
      static constexpr double d1 = -12715105075.0/11282082432.0,  d3 = 87487479700.0/32700410799.0,
                              d4 = -10690763975.0/1880347072.0,   d5 = 701980252875.0/199316789632.0,
                              d6 = -1453857185.0/822651844.0,     d7 = 69997945.0/29380423.0;

      double t_ = 0.0, t_prev_ = 0.0, h_ = 0.0;
      Vec y_, y_new_, y_prev_, tmp_, err_;
      Vec k1_, k2_, k3_, k4_, k5_, k6_, k7_;
      Vec r2_, r3_, r4_, r5_;   // Dense output coefficients of the last accepted step
};


// ====== Driver ======

namespace detail {
   // An event fires when g leaves its sign. Starting exactly at 0 (e.g. right after a restart) does not count.
   inline bool crossed(const double g_prev, const double g) {
      return (g_prev < 0.0 && g >= 0.0) || (g_prev > 0.0 && g <= 0.0);
   }

   // Index of the first component of g that crossed relative to g_prev, -1 if none
   inline int first_crossing(CRef<Vec> g_prev, CRef<Vec> g) {
      for (Eigen::Index i = 0; i < g.size(); ++i)
         if (crossed(g_prev[i], g[i])) return static_cast<int>(i);
      return -1;
   }

   // Bisection on the dense output for the earliest crossing in the last step.
   // On entry g holds the event values at stepper.t() and has a crossing.
   // Returns the time just after the root, so g is already on the new side there.
   template <class Stepper, class Events>
   double locate_event(const Stepper& stepper, Events& events, CRef<Vec> g_prev,
                       Ref<Vec> g, Ref<Vec> g_mid, Ref<Vec> z_mid) {
      double a = stepper.t_prev(), b = stepper.t();   // Invariant: no crossing at a, crossing at b
      const double tol = 4.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(b));

      while (b - a > tol) {
         const double m = 0.5 * (a + b);
         if (m <= a || m >= b) break;

         stepper.interpolate(m, z_mid);
         events(m, z_mid, g_mid);
         if (first_crossing(g_prev, g_mid) >= 0) { b = m; g = g_mid; }
         else a = m;
      }
      return b;
   }

   inline double min_step(const double t) {
      return 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(t));
   }
}

// Integrates z from t0 to tend with any stepper that provides
//   init(rhs, t0, z0), try_step(rhs, h, rtol, atol) -> err, accept(), interpolate(t, out),
//   t(), t_prev(), state() and a static default_controller().
//
// rhs:       void(double t, CRef<Vec> z, Ref<Vec> dz)
// events:    void(double t, CRef<Vec> z, Ref<Vec> g)   g has n_events components
// on_output: bool(double t, CRef<Vec> z)               at t_k = k * dt_out in (t0, stop]; false => stop
// on_step:   bool(double t, CRef<Vec> z)               after every accepted step;        false => stop
//
// Stops at tend, at the first event, or when an observer returns false.
// z is overwritten with the state at the returned time.
template <class Stepper, class Rhs, class Events, class OnOutput, class OnStep>
RKResult integrate(Stepper& stepper, Rhs&& rhs, Events&& events, const int n_events,
                   OnOutput&& on_output, OnStep&& on_step,
                   const double t0, const double tend, Ref<Vec> z, const RKOptions& opt)
{
   PIController controller = Stepper::default_controller();
   stepper.init(rhs, t0, z);

   Vec g_prev(n_events), g(n_events), g_mid(n_events);
   Vec z_stop(z.size()), z_out(z.size());
   events(t0, z, g_prev);

   const bool output = opt.dt_out > 0.0;
   long k_out = output ? static_cast<long>(std::floor(t0 / opt.dt_out)) + 1 : 0;

   double h = std::min({opt.h0, opt.hmax, tend - t0});

   for (int step = 0; step < opt.max_steps; ++step) {
      const double t = stepper.t();
      if (tend - t <= detail::min_step(t)) {
         z = stepper.state();
         return {StopReason::ReachedEnd, t, -1};
      }

      h = std::min(h, tend - t);
      if (h < detail::min_step(t)) {
         z = stepper.state();
         return {StopReason::StepTooSmall, t, -1};
      }

      const double err = stepper.try_step(rhs, h, opt.rtol, opt.atol);
      if (!(err <= 1.0)) {   // Also rejects NaN
         h = controller.reject(err, h);
         continue;
      }
      const double h_next = controller.accept(err, h);
      stepper.accept();

      // Events: the step ends at the earliest root, the rest of the step never happens
      events(stepper.t(), stepper.state(), g);
      int event = detail::first_crossing(g_prev, g);
      double t_stop = stepper.t();
      if (event >= 0) {
         t_stop = detail::locate_event(stepper, events, g_prev, g, g_mid, z_stop);
         event = detail::first_crossing(g_prev, g);
         stepper.interpolate(t_stop, z_stop);
      } else {
         z_stop = stepper.state();
      }

      // Output grid points in (t_prev, t_stop]
      if (output) {
         for (double t_out = k_out * opt.dt_out; t_out <= t_stop; t_out = ++k_out * opt.dt_out) {
            stepper.interpolate(t_out, z_out);
            if (!on_output(t_out, z_out)) {
               z = z_out;
               return {StopReason::ObserverStop, t_out, -1};
            }
         }
      }

      if (!on_step(t_stop, z_stop)) {
         z = z_stop;
         return {StopReason::ObserverStop, t_stop, -1};
      }

      if (event >= 0) {
         z = z_stop;
         return {StopReason::Event, t_stop, event};
      }

      g_prev.swap(g);
      h = std::min(h_next, opt.hmax);
   }

   z = stepper.state();
   return {StopReason::MaxSteps, stepper.t(), -1};
}

// Dormand-Prince 5(4) integration, see integrate()
template <class Rhs, class Events, class OnOutput, class OnStep>
RKResult integrate_rk45(Rhs&& rhs, Events&& events, const int n_events, OnOutput&& on_output, OnStep&& on_step,
              const double t0, const double tend, Ref<Vec> z, const RKOptions& opt)
{
   DP5 stepper(z.size());
   return integrate(stepper, rhs, events, n_events, on_output, on_step, t0, tend, z, opt);
}

} // Namespace rk45
