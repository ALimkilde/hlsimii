#pragma once

#include "eigen_include.h"
#include "line.h"
#include "mesh.h"
#include "linalg.h"

struct Params {
   Vec2 gravity{0.0, -9.81};

   double newton_tol = 1e-5;
   // Drag and dampening parameters
};

struct workspace {
   Vec state; // [x1 y1 x2 y2 ...| vx1 vy1 vx2 vy2 ...]
};

struct PointMass {
   int node;
   double mass;

   PointMass(int n, double m) : node(n), mass(m) {}
};

class LineModel {
  public:
      explicit LineModel(const DiscreteLine& line, 
                         const Vec2 anchor_left, 
                         const Vec2 anchor_right, 
                         const Params& p);  

      Eigen::Index size() const;

      // === Dynamic Solves === //
      void rhs(double t, const Vec& y, Vec& dy) const; // ODE: y input state; dy output rhs = [v | F/m]
                                
      // === Static Solves === //
      void residual(const Vec& q, Vec& r, SymmBlockTriMat* jac = nullptr) const;
      bool static_solver(Vec& q, double tol) const;
      Vec static_solver_initial_guess() const;

      bool both_neighboring_edges_slack_flat(const Vec& q) const;
      bool both_neighboring_edges_slack(Eigen::Ref<const Mat2X> q) const;

      // === Tools to place slackliner ===
      int nearest_node(const Vec& q, const double x_coor) const;
      void place_slackliner(const int node, double mass) const;
      void remove_slackliner() const;

      // === Edge forces === //
      Vec2 edge_vector(Eigen::Ref<const Mat2X> q, Eigen::Index e) const;
      Vec2 edge_force(Eigen::Index e, const Vec2& dq, Mat2* K = nullptr) const; 

      void assemble(Eigen::Ref<const Mat2X> q, Eigen::Ref<Mat2X> F, SymmBlockTriMat* K = nullptr) const; // Compute net forces in F and Jac in K

  private:
      const DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;

      // States related to slackliner
      mutable bool has_slackliner_;
      mutable std::vector<double> m_;        // node masses
      mutable std::vector<double> inv_m_;    // inverse node masses
      mutable int node_slackliner_;
      
};
