#include "physics.h"

LineModel::LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p)
   : line_(std::move(line)),
     anchor_left_(anchor_left),
     anchor_right_(anchor_right),
     params_(p),
     Fe_(2, line.num_elements())
{
}

// Fills edge forces into Fe_
template <bool WithVelocity> // Only add KV dampening if velocity is included
const Vec2 LineModel::edge_force(Eigen::Index e, const Vec2& qa, const Vec2& qb)
{
      Vec2 dq = qb - qa;
      double len = dq.norm();  // Length of e
                               //
      double stretch_main = len - line_.l_main[e];
      double stretch_backup = len - line_.l_backup[e];

      if (stretch_main > 0.0 || stretch_backup > 0.0){
         Vec2 unit = dq/len; // Unit vector

         return unit * ( line_.k_main[e] * std::max(stretch_main, 0.0)
                        + line_.k_backup[e] * std::max(stretch_backup, 0.0) );
      }
      else {
         return Vec2::Zero();
      }

}

// Fills edge forces into Fe_
template <bool WithVelocity> // Only add KV dampening if velocity is included
void LineModel::edge_forces(Eigen::Ref<const Mat2X> q, Eigen::Ref<const Mat2X> v)
{
   const Eigen::Index num_edges = Fe_.cols();

   // Edge e connects global nodes e and e+1; free node k is q.col(k-1).
   Fe_.col(0) = edge_force<WithVelocity>(0, anchor_left_, q.col(0));

   for (Eigen::Index e = 1; e < num_edges - 1; ++e)
      Fe_.col(e) = edge_force<WithVelocity>(e, q.col(e - 1), q.col(e));

   Fe_.col(num_edges - 1) = edge_force<WithVelocity>(num_edges - 1, q.col(num_edges - 2), anchor_right_);

}
