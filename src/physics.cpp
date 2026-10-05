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

// dq = qb - qa for edge e, which runs from node e-1 to node e.
// Node -1 is anchor_left_ and node q.cols() is anchor_right_.
Vec2 LineModel::edge_vector(Eigen::Ref<const Mat2X> q, Eigen::Index e) const {

   const Eigen::Index num_nodes = q.cols();
   assert(e >= 0 && e <= num_nodes);

   const Vec2 qa = (e == 0)         ? anchor_left_  : Vec2(q.col(e - 1));
   const Vec2 qb = (e == num_nodes) ? anchor_right_ : Vec2(q.col(e));

   return qb - qa;
}

  Vec2 LineModel::edge_force(Eigen::Index e, const Vec2& dq, Mat2* K) const
{
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
void LineModel::assemble(Eigen::Ref<const Mat2X> q, Eigen::Ref<Mat2X> F, SymmBlockTriMat* K) const
{
   const Eigen::Index num_edges = line_.num_elements();
   const Eigen::Index num_nodes = q.cols();

   assert(num_edges > 1);

   Mat2 last_edge_K, next_edge_K;
   Mat2* last_edge_Kptr = nullptr; 
   Mat2* next_edge_Kptr = nullptr;
   
   if (K) {
      last_edge_Kptr = &last_edge_K;
      next_edge_Kptr = &next_edge_K;
   }

   Vec2 last_edge_force = edge_force(last_edge(0), edge_vector(q,0), last_edge_Kptr);   
   Vec2 next_edge_force;

   // Todo; refactor an assemble over edges

   for (Eigen::Index i = 0; i < num_nodes; ++i){
      next_edge_force = edge_force(next_edge(i), edge_vector(q,next_edge(i)), next_edge_Kptr);
      F.col(i) = next_edge_force - last_edge_force 
                  + params_.gravity * line_.node_mass[mass_index(i)];

      last_edge_force = next_edge_force;

      if (K) {
         K->diag(i) = -next_edge_K - last_edge_K;
         if (i < num_nodes - 1) K->upper(i) = next_edge_K;

         last_edge_K = next_edge_K;
      }

   }

}

void LineModel::residual(const Vec& q, Vec& r, SymmBlockTriMat* jac) const {

   assert(r.size() == q.size());

   Eigen::Map<const Mat2X> Q(q.data(), 2, q.size()/2);
   Eigen::Map<Mat2X> R(r.data(), 2, r.size()/2);

   assemble(Q, R, jac);

}

Vec static_solver_initial_guess() const {

   // For claude to write
   //
}


bool LineModel::any_edge_slack_flat(const Vec& q) const {

   Eigen::Map<const Mat2X> Q(q.data(), 2, q.size()/2);
   return any_edge_slack(Q);

}


bool LineModel::any_edge_slack(Eigen::Ref<const Mat2X> q) const {

   const Eigen::Index num_edges = line_.num_elements();

   for (Eigen::Index e = 0; e < num_edges; e++) {
      Vec2 dq = edge_vector(q,e);
      const double len = dq.norm();
      const double s_main = len - line_.l_main[e];
      const double s_backup = len - line_.l_backup[e];

      if (s_main<=0 && s_backup<=0) return true;

   }

   return false;

}

void LineModel::static_solver(Vec& q, double tol) const {

   int maxsteps = 1000;
   double reduce_alpha = 0.9;

   Vec r(q.size());
   Vec dq(q.size()), qnew(q.size()); 
   SymmBlockTriMat Jac(q.size()/2);

   for (int i = 0; i < maxsteps; i++) {
       residual(q, r, &Jac);

       if (r.norm() < tol) return

       Jac.solve(dq, -r);

       double alpha = 1;
       qnew = q + alpha*dq;

       // Linesearch to avoid slack edges.
       while ( any_edge_slack_flat(qnew) ){
          alpha *= reduce_alpha;

          qnew = q + alpha*dq;
       }

       q = qnew;

   }

}


