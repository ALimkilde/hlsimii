#pragma once

#include <Eigen/Dense>

#include <line.h>
#include <mesh.h>

using Vec2 = Eigen::Vector2d;
using Vec = Eigen::VectorXd;
using Mat2X = Eigen::Matrix2Xd;

struct Params {
   Vec2 gravity{0.0, -9.81};
   // Drag and dampening parameters
};

struct workspace {
   Vec state; // [x1 y1 x2 y2 ...| vx1 vy1 vx2 vy2 ...]
};

class LineModel {
  public:
      explicit LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p);  

      Eigen::Index size() const;

      // Net node forces F (2 x free nodes) at positions q and velocities v. Fills Fe_ as a side effect.
      void forces(const Eigen::Ref<const Mat2X>& q, const Eigen::Ref<const Mat2X>& v, Mat2X& F);
      // Same as forces() with v = 0: velocity terms (damping, drag) are skipped.
      void static_forces(const Eigen::Ref<const Mat2X>& q, Mat2X& F);

      void rhs(double t, const Vec& y, Vec& dy); // ODE: y input state; dy output rhs = [v | F/m]
      void residual(const Vec& q, Vec& r);       // static equilibrium: r = F(q, 0), in force units

  private:
      template <bool WithVelocity> // Only add KV dampening if velocity is included
      const Vec2 edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb); // Fills Fe_ 
                                                                       //
      template <bool WithVelocity> // Only add KV dampening if velocity is included
      void edge_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v); // Fills Fe_
                                                                                     //
      template <bool WithVelocity> // Only add drag if velocity is included
      void net_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v); // Fills Fe_

      DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;
      mutable Mat2X Fe_;     // per-edge forces 
      Vec inv_m_;    // inverse node masses
};
