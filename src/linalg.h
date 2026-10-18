#pragma once

#include "eigen_include.h"
#include <vector>


// Implement the linear system of Rosenbrock type ODE for general dense Matrices
//     A x = (I/(gamma * h) - J) x = b
class RodasDenseLinSys {
   public:
      explicit RodasDenseLinSys(const Eigen::Index n)
         : n_(n), jac_(Mat::Zero(n,n)), linsys_(n,n), lu_(linsys_) {};

      // return reference for jacobian to be overridden
      Mat& jacobian() { return jac_; }

      // Factorizes
      // (I * fac - J)
      // Eg. in Rodas; call as factorize(1/(gamma*h))
      void factorize(const double fac) 
      { 
         linsys_ = - jac_;
         linsys_.diagonal().array() += fac;

         lu_.compute(linsys_); // Factorize happens here. Overrides linsys.
      };

      // Ax = b
      void solve(Ref<Vec> x, CRef<Vec> b) const 
      {
         x = lu_.solve(b);
      }

      // Getters
      Eigen::Index size() const { return n_; }

   private:
      const Eigen::Index n_;
      Mat jac_, linsys_;
      Eigen::PartialPivLU< Eigen::Ref<Mat> > lu_;
};

// Implements 2x2 block tridiagonal matrix
class BlockTriMat {
   public: 
      // Create zeroed matrix
      explicit BlockTriMat(const int num_blocks);

      void solve(Eigen::Ref<Vec> x, Eigen::Ref<const Vec> b) const;

      void setZero();

      Mat to_dense() const;

      Mat2& upper(const int i) { return U[i]; }
      const Mat2& upper(const int i) const { return U[i]; }

      Mat2& lower(const int i) { return L[i]; }
      const Mat2& lower(const int i) const { return L[i]; }

      Mat2& diag(const int i) { return D[i]; }
      const Mat2& diag(const int i) const { return D[i]; }

      void multiply_rowwise(CRef<Vec> fac);

      int n() const { return num_blocks_ * 2; }
      int num_blocks() const { return num_blocks_; }
   private: 
      int num_blocks_; 
      std::vector<Mat2> D, U, L; // D : Diagonal blocks 
                              // U : is the upper diagonal block
                              // L : is the lower diagonal blocks

};

struct Jac { 
   BlockTriMat& Kq; 
   BlockTriMat& Kv;
};

// Implements the linear system of Rosenbrock type ODE solvers for a highline simulation
//     A x = (I/(gamma * h) - J) x = b
// It considers the special case where x = [q | v] is the reduction 
// of a second order ODE to a first order. This gives the structure
//    A = [ I/gamma*h | -I               ]
//        [ -J_vq     | I/gamma*h - J_vv ]
// Furthermore we assume that J_vq and J_vv is of type BlocKTriMat
class RodasLineModelLinSys {
   public:
      explicit RodasLineModelLinSys(const Eigen::Index n)
         : n_(n), linsys_(n/2,n/2), lu_(linsys_), Kq_(n/4), Kv_(n/4) 
      {
      };

      // return reference for jacobian to be overridden
      BlockTriMat& Kq() { return Kq_; }
      BlockTriMat& Kv() { return Kv_; }
      Jac jacobian() { return {Kq_, Kv_}; }

      // Factorizes
      // (I * fac - J)
      // Eg. in Rodas; call as factorize(1/(gamma*h))
      void factorize(const double fac) 
      { 
         fac_ = fac;

         // TODO; not very efficient
         linsys_ = - Kv_.to_dense() - Kq_.to_dense()/fac;
         linsys_.diagonal().array() += fac;

         lu_.compute(linsys_); // Factorize happens here. Overrides linsys.
      };

      // Ax = b
      void solve(Ref<Vec> x, CRef<Vec> b) const 
      {
         CRef<Vec> bq = b.head(n_/2);
         CRef<Vec> bv = b.tail(n_/2);

         Ref<Vec> xq = x.head(n_/2);
         Ref<Vec> xv = x.tail(n_/2);

         xv = lu_.solve(bv + Kq_.to_dense()*bq/fac_);
         xq = (bq + xv)/fac_;
      }

      // Getters
      Eigen::Index size() const { return n_; }

   private:
      const Eigen::Index n_;
      Mat linsys_;
      Eigen::PartialPivLU< Eigen::Ref<Mat> > lu_;

      double fac_;

      BlockTriMat Kq_, Kv_;
};
