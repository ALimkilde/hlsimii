#include "physics.h"
#include <algorithm>

LineModel::LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p)
   : line_(std::move(line)),
     anchor_left_(anchor_left),
     anchor_right_(anchor_right),
     params_(p),
     Fe_(2, line.num_elements())
{
}

  Vec2 LineModel::edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb, Mat2* K) const
  {
     const Vec2 dq = qb - qa;
     const double len = dq.norm();
     const double s_main = len - line_.l_main[e];
     const double s_backup = len - line_.l_backup[e];
   
     double T = 0.0, dT = 0.0;  // tension and dT/dlen
     if (s_main > 0.0)   { T += line_.k_main[e]   * s_main;   dT += line_.k_main[e]; }
     if (s_backup > 0.0) { T += line_.k_backup[e] * s_backup; dT += line_.k_backup[e]; }
   
     if (T == 0.0) {                 // slack edge
        if (K) K->setZero();
        return Vec2::Zero();
     } 
   
     const Vec2 u = dq / len; // Unit vector
     if (K) {
        const Mat2 uuT = u * u.transpose();
        *K = dT * uuT + (T / len) * (Mat2::Identity() - uuT);  // axial + geometric stiffness
     }  
     return T * u;
  }

// Fills edge forces into Fe_
void LineModel::edge_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v) const
{
   const Eigen::Index num_edges = Fe_.cols();

   Fe_.col(0) = edge_force(0, anchor_left_, q.col(0));

   for (Eigen::Index e = 1; e < num_edges - 1; ++e)
      Fe_.col(e) = edge_force(e, q.col(e - 1), q.col(e));

   Fe_.col(num_edges - 1) = edge_force(num_edges - 1, q.col(num_edges - 2), anchor_right_);

}


