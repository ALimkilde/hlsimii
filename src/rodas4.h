#pragma once

#include "eigen_include.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

// Rodas4: 4th order Rosenbrock method with a 3rd order embedded error estimate.
//         It is based off of campter IV.7 of the book
//
//             "Solving Ordinary Differential Equations II" by E. Hairer & G. Wanner
//
//         and their fortran code "rodas.f" and "dc_decsol.f" from
//
//             https://www.unige.ch/~hairer/software.html
//

namespace rodas4 {

struct Rodas4Options {
   double rtol = 1e-5, atol = 1e-5;
   int max_steps = 1'000'000;
   double h0 = 1e-4, hmax = std::numeric_limits<double>::infinity();
   double dt_out = 0.0;   // Output on the global grid t_k = k * dt_out; <= 0 disables output
   double hmin = 1e-8;
   double time_slack = 1e-8; // Must be >= hmin!
};

enum class StopReason { ReachedEnd, Event, ObserverStop, StepTooSmall, MaxSteps };

struct Rodas4Result {
   StopReason reason;
   double t;          // where integration stopped
   int event_index;   // which event fired (-1 if none)
   int n_accepts;        
   int n_rejects;        
};

// Returns new stepsize
// Assumes err is scaled as err = ||(embedded - ynew)/(atol + rtol |ynew|)||
inline double newStepSize(const double h, double err, bool after_reject = false) {
   const double max_incr = after_reject? 1.0 : 6.0;
   static constexpr double min_incr = 0.2;
   static constexpr double safety   = 0.9;

   if (err == 0.0) return h * max_incr;
   if (!std::isfinite(err)) return h * min_incr;

   const double tolerr   = std::pow(1.0 / err, 0.25); 
   return h * std::min(max_incr, std::max(min_incr, safety * tolerr));
}
   

template <class LinSys>
class Rodas4 {
   public:
      explicit Rodas4(const Eigen::Index n)
         : yprev_(n), y_(n), ynew_(n), dy_(n), tmp_rhs_(n), embedded_(n),
           yinterp1_(n), yinterp2_(n),
           linsys_(n),
           u1_(n), u2_(n), u3_(n), u4_(n), u5_(n), u6_(n) {}

      void init(const double t0, CRef<Vec> z0) {
         t_ = t0;
         tprev_ = t0;
         yprev_ = z0;
         y_ = z0;
         ynew_ = z0;
         jac_is_current_ = false;
      }

      // Autonomous RHS only!
      template <class Rhs, class RhsJac>
      double try_step(Rhs& rhs, RhsJac& rhsjac, const double h, const double rtol, const double atol) {
         h_ = h;

         if (!jac_is_current_) {
            rhsjac(t_, y_, dy_, linsys_.jacobian()); // Fill jac into linsys
            jac_is_current_ = true;
         }

         // Factorize - must not modify jacobian!
         linsys_.factorize(1.0/(gamma * h));
         // ----------------------------------------------------------------

         // Stages
         linsys_.solve(u1_, dy_);

         ynew_ = y_ + a21*u1_;
         rhs(t_ + c2*h, ynew_, tmp_rhs_);
         // TODO; store tmp_rhs + .. in buffers to avoid heap allocations!
         linsys_.solve(u2_, tmp_rhs_ + c21/h * u1_); 

         ynew_ = y_ + a31*u1_ + a32*u2_;
         rhs(t_ + c3*h, ynew_, tmp_rhs_);
         linsys_.solve(u3_, tmp_rhs_ + c31/h * u1_ + c32/h * u2_);

         ynew_ = y_ + a41*u1_ + a42*u2_ + a43*u3_;
         rhs(t_ + c4*h, ynew_, tmp_rhs_);
         linsys_.solve(u4_, tmp_rhs_ + c41/h * u1_ + c42/h * u2_ + c43/h * u3_);

         ynew_ = y_ + a51*u1_ + a52*u2_ + a53*u3_ + a54*u4_;
         rhs(t_ + h, ynew_, tmp_rhs_);
         linsys_.solve(u5_, tmp_rhs_ + c51/h * u1_ + c52/h * u2_ + c53/h * u3_+ c54/h * u4_);

         // Embedded solution
         embedded_ = ynew_ + u5_; // Holds Embedded solution
         rhs(t_ + h, embedded_, tmp_rhs_);
         linsys_.solve(u6_, tmp_rhs_ + c61/h * u1_ + c62/h * u2_ + c63/h * u3_+ c64/h * u4_ + c65/h * u5_);

         ynew_ = embedded_ + u6_; // Finally construct ynew
         
         // normalize error
         u6_.array() /= atol + rtol * y_.array().abs().max(ynew_.array().abs());

         return std::sqrt(u6_.squaredNorm() / static_cast<double>(u6_.size()));
      }

      void accept() {
         yprev_ = y_;
         y_ = ynew_;

         // For dense output
         yinterp1_ = d21*u1_ + d22*u2_ + d23*u3_ + d24*u4_ + d25*u5_; 
         yinterp2_ = d31*u1_ + d32*u2_ + d33*u3_ + d34*u4_ + d35*u5_; 

         tprev_ = t_;
         t_ += h_;

         jac_is_current_ = false;
      }

      void interpolate(const double t, Ref<Vec> out) const {
         static constexpr double hair = 1e-10;
         assert(t >= tprev_ - hair);
         assert(t <= t_ + hair);

         if (t_ == tprev_) {
            out = yprev_; 
            return;
         }

         const double s  = (t - tprev_) / (t_ - tprev_);
         const double s1 = 1.0 - s;

         out = yprev_*s1 + s*(y_ + s1*(yinterp1_ + s*yinterp2_));

      }

      // Getters
      double t() const { return t_; }
      double tprev() const { return tprev_; }
      const Vec& state() const { return y_; } 
      bool jac_is_current() const { return jac_is_current_; }

  private:
      double t_ = 0.0, tprev_ = 0.0, h_ = 0.0;
      Vec yprev_, y_, ynew_, dy_, tmp_rhs_, embedded_; // TODO; get rid of ynew? and more?
      Vec yinterp1_, yinterp2_;
      LinSys linsys_;

      bool jac_is_current_ = false;

      // stages
      Vec u1_, u2_, u3_, u4_, u5_, u6_;

      // Tableau
      static constexpr double c2=0.386e0;
      static constexpr double c3=0.21e0 ;
      static constexpr double c4=0.63e0;
      static constexpr double bet2p=0.0317e0;
      static constexpr double bet3p=0.0635e0;
      static constexpr double bet4p=0.3438e0 ;
      static constexpr double d1=  0.2500000000000000e+00;
      static constexpr double d2= -0.1043000000000000e+00;
      static constexpr double d3=  0.1035000000000000e+00;
      static constexpr double d4= -0.3620000000000023e-01;
      static constexpr double a21= 0.1544000000000000e+01;
      static constexpr double a31= 0.9466785280815826e+00;
      static constexpr double a32= 0.2557011698983284e+00;
      static constexpr double a41= 0.3314825187068521e+01;
      static constexpr double a42= 0.2896124015972201e+01;
      static constexpr double a43= 0.9986419139977817e+00;
      static constexpr double a51= 0.1221224509226641e+01;
      static constexpr double a52= 0.6019134481288629e+01;
      static constexpr double a53= 0.1253708332932087e+02;
      static constexpr double a54=-0.6878860361058950e+00;
      static constexpr double c21=-0.5668800000000000e+01;
      static constexpr double c31=-0.2430093356833875e+01;
      static constexpr double c32=-0.2063599157091915e+00;
      static constexpr double c41=-0.1073529058151375e+00;
      static constexpr double c42=-0.9594562251023355e+01;
      static constexpr double c43=-0.2047028614809616e+02;
      static constexpr double c51= 0.7496443313967647e+01;
      static constexpr double c52=-0.1024680431464352e+02;
      static constexpr double c53=-0.3399990352819905e+02;
      static constexpr double c54= 0.1170890893206160e+02;
      static constexpr double c61= 0.8083246795921522e+01;
      static constexpr double c62=-0.7981132988064893e+01;
      static constexpr double c63=-0.3152159432874371e+02;
      static constexpr double c64= 0.1631930543123136e+02;
      static constexpr double c65=-0.6058818238834054e+01;
      static constexpr double gamma= 0.2500000000000000e+00;

      static constexpr double d21= 0.1012623508344586e+02;
      static constexpr double d22=-0.7487995877610167e+01;
      static constexpr double d23=-0.3480091861555747e+02;
      static constexpr double d24=-0.7992771707568823e+01;
      static constexpr double d25= 0.1025137723295662e+01;
      static constexpr double d31=-0.6762803392801253e+00;
      static constexpr double d32= 0.6087714651680015e+01;
      static constexpr double d33= 0.1643084320892478e+02;
      static constexpr double d34= 0.2476722511418386e+02;
      static constexpr double d35=-0.6594389125716872e+01;
   
};

class OutputGrid {
   public:
      explicit OutputGrid(double dt, double t0) : dt_(dt) {
         disabled_ = (dt <= 0);
         if (disabled_) return;

         i_ = std::floor(t0/dt) + 1;
      }

      double next() const { 
         return disabled_ ? std::numeric_limits<double>::infinity() : i_*dt_; 
      }

      void advance() { i_++; }

   private:
      double dt_;
      Eigen::Index i_; // Where on grid we are
      bool disabled_;

};

// Event handling helper routines
//    - Finding first crossings
//    - Bysection root finding on dense ouput (using stepper.interpolate)
//
// Implemented by Claude
namespace eventHandling {
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
      double a = stepper.tprev(), b = stepper.t();   // Invariant: no crossing at a, crossing at b
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
}

template <class Stepper, class Rhs, class RhsJac, class Events, class OnOutput, class OnStep>
Rodas4Result integrate(Stepper& stepper, Rhs&& rhs, RhsJac rhsjac, Events&& events, const int n_events, OnOutput&& on_output, OnStep&& on_step, const double t0, const double tend, Ref<Vec> z, const Rodas4Options& opt){

   stepper.init(t0, z);

   double h = std::min({opt.h0, opt.hmax, tend - t0}); // step size

   // To store z after stopped event
   double t_stop;
   Vec z_stop(z.size()); 
   
   bool after_reject = false;

   // Event handling
   int event = -1; // -1 means no event happened
   Vec g_prev(n_events), g(n_events), g_mid(n_events);
   events(t0, z, g_prev);

   // To store dense outputs
   OutputGrid grid(opt.dt_out, t0);
   double t_out;
   Vec z_out(z.size());

   int n_rejects = 0;
   int n_accepts = 0;
   
   for (int step = 0; step < opt.max_steps; ++step) {
      const double t = stepper.t();
      if (tend - t <= opt.time_slack) {
         z = stepper.state();

         // Check for one last output
         t_out = tend;
         if (grid.next() <= t_out + opt.time_slack) { // Add a little slack here
            z_out = z;
            if (!on_output(t_out, z_out)) { // Calls on_output
               z = z_out;
               return {StopReason::ObserverStop, t_out, -1, n_accepts, n_rejects};
            }
            grid.advance();
         }

         return {StopReason::ReachedEnd, t, -1, n_accepts, n_rejects};
      }

      h = std::min(h, tend - t);
      if (h < opt.hmin) {
         z = stepper.state();
         return {StopReason::StepTooSmall, t, -1, n_accepts, n_rejects};
      }

      const double err = stepper.try_step(rhs, rhsjac, h, opt.rtol, opt.atol);
      if (!(err <= 1.0)) {   // Also rejects NaN
         h = newStepSize(h, err, after_reject); 
         after_reject = true;
         n_rejects++; 
         continue;
      }
      
      const double h_next = newStepSize(h, err, after_reject); 
      stepper.accept();
      after_reject = false;
      n_accepts++;

      // Event handling
      events(stepper.t(), stepper.state(), g);
      event = eventHandling::first_crossing(g_prev, g);
      if (event >= 0) {
         t_stop = eventHandling::locate_event(stepper, events, g_prev, g, g_mid, z_stop);
         event = eventHandling::first_crossing(g_prev, g);
         stepper.interpolate(t_stop, z_stop);
      } else { // No event
         z_stop = stepper.state();
         t_stop = stepper.t();
      }

      // Dense output
      while(grid.next() <= t_stop) {
         t_out = grid.next();
         stepper.interpolate(t_out, z_out);
         if (!on_output(t_out, z_out)) { // Calls on_output
            z = z_out;
            return {StopReason::ObserverStop, t_out, -1, n_accepts, n_rejects};
         }

         grid.advance();
      }

      if (!on_step(t_stop, z_stop)) { // Note that this calls on_step for every step
         z = z_stop;
         return {StopReason::ObserverStop, t_stop, -1, n_accepts, n_rejects};
      }

      if (event >= 0) {
         z = z_stop;
         return {StopReason::Event, t_stop, event, n_accepts, n_rejects};
      }

      // Update event
      g_prev.swap(g);
      // Step size for next step
      h = std::min(h_next, opt.hmax);
   }

   z = stepper.state();
   return {StopReason::MaxSteps, stepper.t(), -1, n_accepts, n_rejects};

};

// Rosenbrock 4(3) integration
template <class LinSys, class Rhs, class RhsJac, class Events, class OnOutput, class OnStep>
Rodas4Result integrate_rodas4(Rhs&& rhs, RhsJac&& rhsjac, Events&& events, const int n_events, OnOutput&& on_output, OnStep&& on_step,
              const double t0, const double tend, Ref<Vec> z, const Rodas4Options& opt)
{
   Rodas4<LinSys> stepper(z.size());
   return integrate(stepper, rhs, rhsjac, events, n_events, on_output, on_step, t0, tend, z, opt);
}

} // Namespace

 
