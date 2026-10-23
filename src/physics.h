#pragma once

#include "eigen_include.h" // Note that Ref is Eigen::Ref<..> and CRef Eigen::Ref<const ..>
#include "line.h"
#include "mesh.h"
#include "linalg.h"

// Views of the ODE state z = [q | v]. q and v are 2 x n_nodes each.
struct StateView {
   Eigen::Map<Mat2X> q;
   Eigen::Map<Mat2X> v;
};

struct ConstStateView {
   Eigen::Map<const Mat2X> q;
   Eigen::Map<const Mat2X> v;
};

inline StateView state_view(Ref<Vec> z) {
   assert(z.size() % 4 == 0);
   const Eigen::Index half = z.size() / 2;
   return { as_mat2x(z.head(half)), as_mat2x(z.tail(half)) };
}

inline ConstStateView const_state_view(CRef<Vec> z) {
   assert(z.size() % 4 == 0);
   const Eigen::Index half = z.size() / 2;
   return { as_const_mat2x(z.head(half)), as_const_mat2x(z.tail(half)) };
}

namespace constants {
   inline const Vec2 gravity{0.0, -9.81};
}

struct Params {
   Vec2 gravity = constants::gravity;
   double tau_KV = 0.23;
   double c_drag = 1.15; 
   double webbing_width = 0.025; 
   double rho_air = 1.225; // Density of air
   double newton_tol = 1e-5;
   // Drag and dampening parameters
};

// Shared data between local assembly (edge_force, edge_drag)
struct EdgeGeometry {
   double len;
   double leninv; // 1/len
   Vec2 u;
   Vec2 n;

   Mat2 uuT;
   Mat2 IuuT; // I - uuT

   explicit EdgeGeometry(const Vec2& dq, bool with_matricies) : 
      len(dq.norm()), 
      leninv(1/len),
      u(dq*leninv),
      n(-u.y(), u.x())
   {
      if (with_matricies) {
         uuT = u*u.transpose();
         IuuT = Mat2::Identity() - uuT;
      }
   }
};

class LineModel {
  public:
      explicit LineModel(const DiscreteLine& line, 
                         const Vec2 anchor_left, 
                         const Vec2 anchor_right, 
                         const Params& p);  

      // === Dynamic Solves === //
      // ODE: z = [q | v} input state; dz output rhs = [v | F/m]
      void rhs(double t, CRef<Vec> z, Ref<Vec> dz, BlockTriMat* Kq = nullptr, BlockTriMat* Kv = nullptr) const;  // Out facing API
      void rhs(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> dq, Ref<Mat2X> dv, BlockTriMat* Kq = nullptr, BlockTriMat* Kv = nullptr) const;  // Internal implementation
                                
      // === Static Solves === //
      void residual(const Vec& q_flat, Vec& r_flat, BlockTriMat* jac = nullptr) const;
      bool static_solver(Vec& q_flat, double tol) const;
      Vec static_solver_initial_guess() const;

      bool both_neighboring_edges_slack_flat(const Vec& q) const;
      bool both_neighboring_edges_slack(CRef<Mat2X> q) const;

      // === Tools to place slackliner ===
      int nearest_node(const Vec& q_flat, const double x_coor) const;
      void place_slackliner(const int node, double mass) const;
      void remove_slackliner() const;

      // === Edge forces === //
      Vec2 edge_vector(CRef<Mat2X> q, Eigen::Index e) const;
      Vec2 edge_force(Eigen::Index e, const EdgeGeometry& geom, const Vec2& dv, Mat2* Kq = nullptr, Mat2* Kv = nullptr) const; 
      Vec2 edge_force(Eigen::Index e, const Vec2& dq, const Vec2& dv, Mat2* Kq = nullptr, Mat2* Kv = nullptr) const; 

      // === KV dampening === //
      // Implementation of KV is in edge_force
      Vec2 vel_vector(CRef<Mat2X> v, Eigen::Index e) const; 

      // === drag === //
      Vec2 edge_drag(Eigen::Index e, const EdgeGeometry& geom, const Vec2& v, Mat2* Kq = nullptr, Mat2* Kv = nullptr) const;
      Vec2 edge_drag(Eigen::Index e, const Vec2& dq, const Vec2& v, Mat2* Kq = nullptr, Mat2* Kv = nullptr) const;

      void assemble(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> F, BlockTriMat* Kq = nullptr, BlockTriMat* Kv = nullptr) const; // Compute net forces in F and Jac in K
      Vec2 node_force(CRef<Mat2X> q, CRef<Mat2X> v, Eigen::Index i) const; // Same as assemble().col(i)

      double get_mass(const int i) const {
         return m_[i];
      }
 
      double gravitational_energy(CRef<Mat2X> q) const;
      double elastic_energy(CRef<Mat2X> q) const;
      double kinetic_energy(CRef<Mat2X> v) const;

      int num_nodes() const { return line_.num_nodes(); }
      int num_edges() const { return line_.num_nodes() + 1; }


  private:
      const DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;


      // States related to slackliner
      mutable bool has_slackliner_;
      mutable int node_slackliner_;

      mutable Vec m_, minv_;        // node masses and inverses
      std::vector<double> drag_coefs_;

      // Compile-time versions with and without Jacobians (public functions dispatch on nullptr), defined in physics.cpp
      template<bool WithKq, bool WithKv> Vec2 edge_force_impl(Eigen::Index e, const EdgeGeometry& geom, const Vec2& dv, Mat2* Kq, Mat2* Kv) const;
      template<bool WithKq, bool WithKv> Vec2 edge_drag_impl(Eigen::Index e, const EdgeGeometry& geom, const Vec2& v, Mat2* Kq, Mat2* Kv) const;
      template<bool WithKq, bool WithKv> void assemble_impl(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> F, BlockTriMat* Kq, BlockTriMat* Kv) const;

};
