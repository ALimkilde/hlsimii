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

void BlockTriMat::multiply_rowwise(CRef<Vec> fac) {

   for(Eigen::Index b = 0; b < num_blocks_; b++) {
      const Eigen::Index u = b;
      const Eigen::Index l = b-1;
      const bool u_exist = u < num_blocks_ - 1;
      const bool l_exist = l >= 0;

      diag(b) *= fac(b);
      if (u_exist) upper(u) *= fac(b);
      if (l_exist) lower(l) *= fac(b);
   }

}

Vec BlockTriMat::operator*(CRef<Vec> x) const
{
   Vec out = Vec::Zero(x.size());

   for(Eigen::Index b = 0; b < num_blocks(); b++) {
      const Eigen::Index u = b;
      const Eigen::Index l = b-1;
      const bool u_exist = u < num_blocks() - 1;
      const bool l_exist = l >= 0;

      out.segment<2>(2*b) += diag(b)*x.segment<2>(2*b);

      if (u_exist) out.segment<2>(2*b) += upper(u)*x.segment<2>(2*(b+1));
      if (l_exist) out.segment<2>(2*b) += lower(l)*x.segment<2>(2*(b-1));
   }
   
   return out;
}

void BlockTriLU::compute(const BlockTriMat& A)
{
   D[0] = A.diag(0).inverse();
   if (1 < num_blocks_) U[0] = D[0]*A.upper(0);

   for (int b = 1; b < num_blocks_; b++) {
      const Eigen::Index u = b;
      const Eigen::Index l = b-1;
      const bool u_exist = u < num_blocks_ - 1;

      D[b] = (A.diag(b) - A.lower(l)*U[b-1]).inverse();
      L[l] = D[b] * A.lower(l);
      if (u_exist) U[b] = D[b]*A.upper(b);
   }
}

void BlockTriLU::solve(Ref<Vec> x, CRef<Vec> b) const
{
   // Forward pass
   d[0] = D[0]*b.segment<2>(0);

   for (int i = 1; i < num_blocks_; i++) {
      const Eigen::Index l = i-1;

      d[i] = D[i]*b.segment<2>(2*i) - L[l]*d[i-1]; 
   }

   // Back substitution
   x.segment<2>(2*(num_blocks_-1)) = d[num_blocks_-1];

   for (int i = num_blocks_-2; i >= 0; i--) {
      x.segment<2>(2*i) = d[i] - U[i]*x.segment<2>(2*(i+1));
   }
}

static void rodasReducedSystem(const BlockTriMat& Kq, const BlockTriMat& Kv, const double fac, BlockTriMat& out) 
{
   for(Eigen::Index b = 0; b < out.num_blocks(); b++) {
      const Eigen::Index u = b;
      const Eigen::Index l = b-1;
      const bool u_exist = u < out.num_blocks() - 1;
      const bool l_exist = l >= 0;

      out.diag(b) = -Kv.diag(b) - Kq.diag(b)/fac;
      out.diag(b).diagonal().array() += fac;

      if (u_exist) out.upper(u) = -Kv.upper(u) - Kq.upper(u)/fac;
      if (l_exist) out.lower(l) = -Kv.lower(l) - Kq.lower(l)/fac;
   }
};

void RodasLineModelLinSys::factorize(const double fac)
{ 
   fac_ = fac;

   rodasReducedSystem(Kq_, Kv_, fac, reduced_system_);

   reduced_lu_.compute(reduced_system_); // Factorize happens here. Overrides linsys.
};

void RodasLineModelLinSys::solve(Ref<Vec> x, CRef<Vec> b) const 
{
   CRef<Vec> bq = b.head(n_/2);
   CRef<Vec> bv = b.tail(n_/2);

   Ref<Vec> xq = x.head(n_/2);
   Ref<Vec> xv = x.tail(n_/2);

   reduced_lu_.solve(xv, bv + Kq_*bq/fac_);
   xq = (bq + xv)/fac_;
}
