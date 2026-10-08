#pragma once

#include "eigen_include.h" // Note that Ref is Eigen::Ref<..> and CRef Eigen::Ref<const ..>
#include "line.h"
#include "mesh.h"
#include "linalg.h"

namespace constants {
   inline const Vec2 gravity{0.0, -9.81};
}

struct Params {
   Vec2 gravity = constants::gravity;
   double newton_tol = 1e-5;
   // Drag and dampening parameters
};

class LineModel {
  public:
      explicit LineModel(const DiscreteLine& line, 
                         const Vec2 anchor_left, 
                         const Vec2 anchor_right, 
                         const Params& p);  

      // === Dynamic Solves === //
      // ODE: z = [q | v} input state; dz output rhs = [v | F/m]
      void rhs(double t, const Vec& z, Vec& dz) const;  // Out facing API
      void rhs(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> dq, Ref<Mat2X> dv) const;  // Internal implementation
                                
      // === Static Solves === //
      void residual(const Vec& q, Vec& r, SymmBlockTriMat* jac = nullptr) const;
      bool static_solver(Vec& q, double tol) const;
      Vec static_solver_initial_guess() const;

      bool both_neighboring_edges_slack_flat(const Vec& q) const;
      bool both_neighboring_edges_slack(CRef<Mat2X> q) const;

      // === Tools to place slackliner ===
      int nearest_node(const Vec& q, const double x_coor) const;
      void place_slackliner(const int node, double mass) const;
      void remove_slackliner() const;

      // === Edge forces === //
      Vec2 edge_vector(CRef<Mat2X> q, Eigen::Index e) const;
      Vec2 edge_force(Eigen::Index e, const Vec2& dq, Mat2* K = nullptr) const; 

      void assemble(CRef<Mat2X> q, Ref<Mat2X> F, SymmBlockTriMat* K = nullptr) const; // Compute net forces in F and Jac in K

  private:
      const DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;

      // States related to slackliner
      mutable bool has_slackliner_;
      mutable int node_slackliner_;

      mutable std::vector<double> m_;        // node masses
      //TODO save inverses of mass to optimize. 
      //Requires updating for place/removing slackliner
      // mutable std::vector<double> inv_m_;    // inverse node masses
      
      
};
