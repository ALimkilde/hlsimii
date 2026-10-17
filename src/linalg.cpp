#include "linalg.h"

BlockTriMat::BlockTriMat(const int num_blocks)
   : num_blocks_(num_blocks),
     D(num_blocks, Mat2::Zero()),
     U(num_blocks - 1, Mat2::Zero()),
     L(num_blocks - 1, Mat2::Zero())
{  
}   

Mat BlockTriMat::to_dense() const {

   Mat out = Mat::Zero(n(), n());

   for (int i = 0; i < num_blocks_; i++){ //block index
      int j = 2*i; //dense index

      out.block<2,2>(j,j) = diag(i);

      if (i < num_blocks_-1) out.block<2,2>(j,j+2) = upper(i);

      if (i > 0) out.block<2,2>(j,j-2) = lower(i-1);
   }

   return out;
}

// Dense reference solve of A x = b via Eigen's LU. O(n^3); for development only.
void BlockTriMat::solve(Eigen::Ref<Vec> x, Eigen::Ref<const Vec> b) const
{
   x = to_dense().partialPivLu().solve(b);
}

void BlockTriMat::setZero() 
{
   for(Eigen::Index b = 0; b < num_blocks_; b++){
      D[b].setZero();
      if (b < num_blocks_ - 1) {
         U[b].setZero();
         L[b].setZero();
      }
   }

}
