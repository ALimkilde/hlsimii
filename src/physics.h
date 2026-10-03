#pragma once

#include <Eigen/Dense>

#include "line.h"
#include "mesh.h"

using Vec2 = Eigen::Vector2d;
using Vec = Eigen::VectorXd;
using Mat2X = Eigen::Matrix2Xd;
using Mat2 = Eigen::Matrix2d;

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

      // === Dynamic Solves === //
      void rhs(double t, const Vec& y, Vec& dy) const; // ODE: y input state; dy output rhs = [v | F/m]
                                
      // === Static Solves === //
      void residual(const Vec& q, Vec& r) const;       // static equilibrium: r = F(q, 0), in force units

      // === Edge forces === //
      Vec2 edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb, Mat2* K = nullptr) const; // Fills Fe_ 
                                                                       //
      void edge_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v) const; // Fills Fe_
                                                                                     //
      void net_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v, Eigen::Ref<Mat2X> F) const; // return sum of forces
                                                                                   //
      Mat2X get_Fe() const { return Fe_; };

  private:
      const DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;
      mutable Mat2X Fe_;     // per-edge forces 
      Vec inv_m_;    // inverse node masses
};
