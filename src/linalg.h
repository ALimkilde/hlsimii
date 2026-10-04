#pragma once

#include "eigen_include.h"
#include <vector>

// Implements symmetric 2x2 block matrix
class SymmBlockTriMat {
   public: 
      // Create zeroed matrix
      explicit SymmBlockTriMat(const int num_blocks);

      void solve(Eigen::Ref<Vec> x, Eigen::Ref<const Vec> b) const;

      Mat to_dense() const;

      Mat2& upper(const int i) { return U[i]; }
      const Mat2& upper(const int i) const { return U[i]; }

      Mat2& diag(const int i) { return D[i]; }
      const Mat2& diag(const int i) const { return D[i]; }

      int n() const { return num_blocks_ * 2; }
      int num_blocks() const { return num_blocks_; }
   private: 
      int num_blocks_; 
      std::vector<Mat2> D, U; // D : Diagonal blocks 
                              // U : is the upper diagonal block
                              // lower diagonal block: U.transpose() 

};
