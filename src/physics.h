#pragma once

#include "eigen_include.h"
#include "line.h"
#include "mesh.h"
#include "linalg.h"

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
      Vec2 edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb, Mat2* K = nullptr) const; 
                                                                                     //
      void assemble(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v, Eigen::Ref<Mat2X> F, SymmBlockTriMat* K = nullptr) const; // Compute net forces in F and Jac in K

  private:
      const DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;
      Vec inv_m_;    // inverse node masses
};
