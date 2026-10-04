#include "physics.h"
#include <algorithm>
#include <cassert>

LineModel::LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p)
   : line_(std::move(line)),
     anchor_left_(anchor_left),
     anchor_right_(anchor_right),
     params_(p)
{
}

  Vec2 LineModel::edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb, Mat2* K) const
{
   const Vec2 dq = qb - qa;
   const double len = dq.norm();
   const double s_main = len - line_.l_main[e];
   const double s_backup = len - line_.l_backup[e];
 
   double T = 0.0, dT = 0.0;   // tension, and dT/dlen
   if (s_main > 0.0)   { 
      T += line_.k_main[e] * s_main;   
      dT += line_.k_main[e]; 
   } 
   if (s_backup > 0.0) { 
      T += line_.k_backup[e] * s_backup; 
      dT += line_.k_backup[e]; 
   }
     
   // === Jacobian ===
   if (T == 0.0) {                 // slack edge
      if (K) K->setZero();
      return Vec2::Zero();
   } 
 
   const Vec2 u = dq / len; // Unit vector
   if (K) {
      Mat2 uuT = u * u.transpose();
      *K = T/len * (Mat2::Identity() - uuT) + dT * uuT;
   }
   
   // === Return rhs ===
   return T * u;
}

static Eigen::Index mass_index(Eigen::Index i) { return i + 1; }
static Eigen::Index last_edge(Eigen::Index i) { return i; }
static Eigen::Index next_edge(Eigen::Index i) { return i+1; }

// Fills net forces into F (and Jacobian if argument present)
void LineModel::assemble(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v, Eigen::Ref<Mat2X> F, SymmBlockTriMat* K) const
{
   const Eigen::Index num_edges = line_.num_elements();

   assert(num_edges > 1);

   Mat2 last_edge_K, next_edge_K;
   Mat2* last_edge_Kptr = nullptr; 
   Mat2* next_edge_Kptr = nullptr;
   
   if (K) {
      last_edge_Kptr = &last_edge_K;
      next_edge_Kptr = &next_edge_K;
   }

   Vec2 last_edge_force = edge_force(last_edge(0), anchor_left_, q.col(0), last_edge_Kptr);   
   Vec2 next_edge_force;


   for (Eigen::Index i = 0; i < num_edges - 2; ++i){
      next_edge_force = edge_force(next_edge(i), q.col(i), q.col(i+1), next_edge_Kptr);
      F.col(i) = next_edge_force - last_edge_force 
                  + params_.gravity * line_.node_mass[mass_index(i)];

      last_edge_force = next_edge_force;

      if (K) {
         K->diag(i) = -next_edge_K - last_edge_K;
         K->upper(i) = next_edge_K;

         last_edge_K = next_edge_K;
      }

   }

   Eigen::Index i = num_edges - 2;
   next_edge_force = edge_force(next_edge(i), q.col(i), anchor_right_, next_edge_Kptr);
   F.col(i) = next_edge_force - last_edge_force 
               + params_.gravity * line_.node_mass[mass_index(i)];

   if (K) K->diag(i) = -next_edge_K - last_edge_K;

}


