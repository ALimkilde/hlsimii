#include "physics.h"
#include <algorithm>
#include <cassert>
#include <iostream>

LineModel::LineModel(const DiscreteLine& line, 
                     const Vec2 anchor_left, 
                     const Vec2 anchor_right, 
                     const Params& p)
   : line_(std::move(line)),
     anchor_left_(anchor_left),
     anchor_right_(anchor_right),
     params_(p),
     has_slackliner_(false),
     node_slackliner_(-1),
     m_(line_.node_mass)
{
}

// TODO move to driver!
int LineModel::nearest_node(const Vec& q_flat, const double x_coor) const {

   // Nearest free node in x to the slackliner
   Eigen::Index nearest;
   (as_const_mat2x(q_flat).row(0).array() - x_coor).abs().minCoeff(&nearest);
   return static_cast<int>(nearest);

}

void LineModel::place_slackliner(const int node, double mass) const {
   if (has_slackliner_) remove_slackliner();

   assert(node >= 0 && node < line_.num_nodes());
   assert(mass > 0);

   has_slackliner_ = true;
   node_slackliner_ = node;
   m_[node_slackliner_] = line_.node_mass[node_slackliner_] + mass;
}

void LineModel::remove_slackliner() const {
   if (!has_slackliner_) return;

   has_slackliner_ = false;
   m_[node_slackliner_] = line_.node_mass[node_slackliner_];
}

// dq = qb - qa for edge e, which runs from node e-1 to node e.
// Node -1 is anchor_left_ and node q.cols() is anchor_right_.
Vec2 LineModel::edge_vector(CRef<Mat2X> q, Eigen::Index e) const {

   const Eigen::Index num_nodes = q.cols();
   assert(e >= 0 && e <= num_nodes);

   const Vec2 qa = (e == 0)         ? anchor_left_  : Vec2(q.col(e - 1));
   const Vec2 qb = (e == num_nodes) ? anchor_right_ : Vec2(q.col(e));

   return qb - qa;
}

Vec2 LineModel::vel_vector(CRef<Mat2X> v, Eigen::Index e) const {

   const Eigen::Index num_nodes = v.cols();
   assert(e >= 0 && e <= num_nodes);

   const Vec2 va = (e == 0)         ? Vec2::Zero()  : Vec2(v.col(e - 1));
   const Vec2 vb = (e == num_nodes) ? Vec2::Zero() : Vec2(v.col(e));

   return vb - va;
}

  Vec2 LineModel::edge_force(Eigen::Index e, const Vec2& dq, const Vec2& dv, Mat2* Kq, Mat2* Kv) const
{
   const double len = dq.norm();
   const double s_main = len - line_.l_main[e];
   const double s_backup = len - line_.l_backup[e];
 
   double T = 0.0, dT = 0.0;   // tension, and dT/dlen
   double c_KV = 0.0;          // Kelvin Voight Dampening constant
                               //
   // TODO; optimization idea; this check doesn't require square root of dq.norm() unless true
   if (s_main > 0.0)   { 
      T += line_.k_main[e] * s_main;   
      dT += line_.k_main[e]; 
      c_KV += line_.k_main[e]; 
   } 
   if (s_backup > 0.0) { 
      T += line_.k_backup[e] * s_backup; 
      dT += line_.k_backup[e]; 
      c_KV += line_.k_backup[e]; 
   }
   c_KV *= params_.tau_KV; // c_KV = tau_KV * k
     
   if (T == 0.0) {                 // slack edge
      if (Kq) Kq->setZero();
      if (Kv) Kv->setZero();
      return Vec2::Zero();
   } 
 
   const Vec2 u = dq / len; // Unit vector

   // Kelvin Voight Dampening
   double feKV = c_KV * dv.dot(u);
   // Clip:
   if (feKV <= -T) {
      feKV = -T;
      if (Kq) Kq->setZero();
      if (Kv) Kv->setZero();
   } else {
      
      if (feKV >= 0) {
         feKV = 0; // Dont dampen stretching edges
         if (Kv) Kv->setZero();
      }

      Mat2 uuT;
      if (Kq || Kv) uuT = u * u.transpose(); // Both derivatives use this

      if (Kq) {
         Mat2 IuuT = (Mat2::Identity() - uuT);
         *Kq = T/len * IuuT + dT * uuT;
         if (feKV < 0) {
            *Kq += 1/len * c_KV * u * dv.transpose() * IuuT; 
            *Kq += 1/len * feKV * IuuT;

         }
      }

      if (!(feKV>=0) && Kv) *Kv = c_KV * uuT; 

   }


   return (T + feKV) * u;
}

// TODO; optimization possible; let this method reuse u, len, uuT, I-uuT.
Vec2 LineModel::edge_drag(Eigen::Index e, const Vec2& dq, const Vec2& v, Mat2* Kq, Mat2* Kv) const {
   const double len = dq.norm();
   const Vec2 u = dq / len; // Unit vector
   const Vec2 vtilde = v - v.dot(u)*u; // Orthogonal part of v

   const double vnorm2 = vtilde.squaredNorm();
   if (vnorm2 == 0.0){
      if (Kq) Kq->setZero();
      if (Kv) Kv->setZero();
      return Vec2::Zero();
   }

   const double c = 0.5 * params_.rho_air * params_.c_drag * params_.webbing_width * (line_.l_main[e] + line_.l_backup[e]);

   const double vlen = std::sqrt(vnorm2);
   Vec2 drag = - c * vtilde * vlen;

   Mat2 uuT, IuuT, I, ddrag_dvtilde;
   if (Kq || Kv) {
      uuT = u * u.transpose(); // Both derivatives use this
      I = Mat2::Identity();
      IuuT = (I - uuT);

      ddrag_dvtilde = c * (-I*vlen - vtilde*vtilde.transpose()/vlen);
   }

   if (Kq) {
      Mat2 dvtilde_dq = (-u*v.transpose() - v.dot(u)*I)*IuuT/len;
      *Kq = ddrag_dvtilde * dvtilde_dq;
   }

   if (Kv) {
      *Kv = ddrag_dvtilde * IuuT;
   }

   return drag;

}

void LineModel::assemble(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> F, BlockTriMat* Kq, BlockTriMat* Kv) const
{
   const Eigen::Index num_edges = line_.num_elements();
   const Eigen::Index num_nodes = q.cols();
   const Vec2& g = params_.gravity;

   Mat2 edge_Kq, drag_Kq;
   Mat2 *edge_Kqptr = nullptr, *drag_Kqptr = nullptr; 
   if (Kq) {
      edge_Kqptr = &edge_Kq;
      drag_Kqptr = &drag_Kq;
      Kq->setZero();
   }

   Mat2 edge_Kv, drag_Kv;
   Mat2 *edge_Kvptr = nullptr, *drag_Kvptr = nullptr; 
   if (Kv) {
      edge_Kvptr = &edge_Kv;
      drag_Kvptr = &drag_Kv;
      Kv->setZero();
   }

   assert(num_edges > 1);
   assert(num_edges == num_nodes + 1);

   for (Eigen::Index i = 0; i < num_nodes; ++i){
      F.col(i) = m_[i] * g;
   }

   for (Eigen::Index e = 0; e < num_edges; ++e){
      const Eigen::Index a = e - 1;          // left end of edge e
      const Eigen::Index b = e;              // right end of edge e
      const bool a_free = (a >= 0);          // is a a free node
      const bool b_free = (b < num_nodes);   // is b a free node

      const Vec2 dq = edge_vector(q, e);
      const Vec2 f = edge_force(e, dq, vel_vector(v, e), edge_Kqptr, edge_Kvptr);

      if (a_free) {
         F.col(a) += f + 0.5 * edge_drag(e, dq, v.col(a), drag_Kqptr, drag_Kvptr);
         if (Kv) Kv->diag(a) += 0.5 * drag_Kv;
         if (Kq) {
            Kq->diag(a) -= 0.5 * drag_Kq;
            if (b_free) {
               Kq->upper(a) += 0.5 * drag_Kq;
            }
         }
      }
      if (b_free) {
         F.col(b) += 0.5 * edge_drag(e, dq, v.col(b), drag_Kqptr, drag_Kvptr) - f;
         if (Kv) Kv->diag(b) += 0.5 * drag_Kv;

         if (Kq) {
            Kq->diag(b) += 0.5 * drag_Kq;

            if (a_free) {
               Kq->lower(a) -= 0.5 * drag_Kq;
            }
         }
      }

      if (Kq) {
         if (a_free) Kq->diag(a) -= edge_Kq;
         if (b_free) Kq->diag(b) -= edge_Kq;
         if (a_free && b_free) Kq->upper(a) += edge_Kq;
         if (a_free && b_free) Kq->lower(a) += edge_Kq;
      }

      if (Kv) {
         if (a_free) Kv->diag(a) -= edge_Kv;
         if (b_free) Kv->diag(b) -= edge_Kv;
         if (a_free && b_free) Kv->upper(a) += edge_Kv;
         if (a_free && b_free) Kv->lower(a) += edge_Kv;
      }

   }

}

void LineModel::residual(const Vec& q_flat, Vec& r_flat, BlockTriMat* jac) const {

   assert(r_flat.size() == q_flat.size());

   Vec v_flat = Vec::Zero(q_flat.size()); // Zero velocity to avoid dampening

   assemble(as_const_mat2x(q_flat), as_const_mat2x(v_flat), as_mat2x(r_flat), jac);

}

// Generated by Claude:
// 
// A non-slack starting shape for the static solver. Every edge gets length
// stretch * l_e, where l_e is the shorter of its main and backup rest lengths.
// If the span is long enough the nodes lie on the straight chord; otherwise the
// chain forms a V hanging below the chord, with the kink placed exactly on a node
// so every edge lies on a straight leg and keeps its full length.
Vec LineModel::static_solver_initial_guess() const {

   const double stretch = 1.01;
   const Eigen::Index num_edges = line_.num_elements();
   const Eigen::Index num_nodes = num_edges - 1;   // free nodes, anchors excluded

   std::vector<double> l(num_edges);
   for (Eigen::Index e = 0; e < num_edges; e++)
      l[e] = std::min(line_.l_main[e], line_.l_backup[e]);

   const Vec2 chord = anchor_right_ - anchor_left_;
   const double D = chord.norm();
   const Vec2 u = chord / D;

   double total = 0.0;
   for (double le : l) total += le;

   Vec q_flat(2 * num_nodes);
   Eigen::Map<Mat2X> q = as_mat2x(q_flat);

   // Straight chord: stretch the chain uniformly to fit the span
   if (D >= stretch * total) {
      double s = 0.0;
      for (Eigen::Index i = 0; i < num_nodes; i++) {
         s += l[i];
         q.col(i) = anchor_left_ + (s / total) * chord;
      }
      return q_flat;
   }

   // V shape: kink at the node closest to half the rest length
   Eigen::Index kink = 0;
   double A = 0.0;
   while (kink < num_nodes - 1 && A + l[kink] < 0.5 * total) A += l[kink++];
   A = stretch * (A + l[kink]);            // left leg: edges 0..kink
   const double B = stretch * total - A;   // right leg: the remaining edges
   assert(std::abs(A - B) <= D);           // the triangle (A, B, D) must exist

   // Kink position from the triangle: a along the chord, h perpendicular to it
   const double a = (A * A - B * B + D * D) / (2.0 * D);
   const double h = std::sqrt(A * A - a * a);
   Vec2 n(-u.y(), u.x());
   if (n.dot(params_.gravity) < 0.0) n = -n;   // hang towards gravity
   const Vec2 apex = anchor_left_ + a * u + h * n;

   double s = 0.0;
   for (Eigen::Index i = 0; i < num_nodes; i++) {
      s += stretch * l[i];
      q.col(i) = (i <= kink) ? Vec2(anchor_left_ + (s / A) * (apex - anchor_left_))
                             : Vec2(apex + ((s - A) / B) * (anchor_right_ - apex));
   }
   return q_flat;
}


bool LineModel::both_neighboring_edges_slack_flat(const Vec& q_flat) const {

   return both_neighboring_edges_slack(as_const_mat2x(q_flat));

}


bool LineModel::both_neighboring_edges_slack(CRef<Mat2X> q) const {

   const Eigen::Index num_edges = line_.num_elements();
   bool last_edge_slack = false;

   for (Eigen::Index e = 0; e < num_edges; e++) {
      // TODO; switch this for calling edge_tension(e, len)
      double T = edge_force(e, edge_vector(q,e), Vec2::Zero()).norm();

      bool this_edge_slack = T <= 1E-12;

      if (last_edge_slack && this_edge_slack) return true;

      last_edge_slack = this_edge_slack;
   }

   return false;

}

bool LineModel::static_solver(Vec& q, double tol) const {

   int maxsteps = 1000;
   double reduce_alpha = 0.9;

   Vec r(q.size());
   Vec dq(q.size()), qnew(q.size()); 
   BlockTriMat Jac(q.size()/2);

   for (int i = 0; i < maxsteps; i++) {
       residual(q, r, &Jac);

       if (r.norm() < tol) { 
          return true; }

       Jac.solve(dq, -r);

       double alpha = 1;
       qnew = q + alpha*dq;

       // Linesearch to avoid slack edges.
       while ( both_neighboring_edges_slack_flat(qnew) ){
          alpha *= reduce_alpha;

          if (alpha < 1e-10) return false; // The solver didn't take a step

          qnew = q + alpha*dq;
       }

       q = qnew;

   }

   return false;

}


void LineModel::rhs(double t, CRef<Vec> z, Ref<Vec> dz, BlockTriMat* Kq, BlockTriMat* Kv) const {

   assert(z.size() ==  dz.size());

   ConstStateView s = const_state_view(z);
   StateView ds = state_view(dz);

   rhs(s.q, s.v, ds.q, ds.v, Kq, Kv);
}


void LineModel::rhs(CRef<Mat2X> q, CRef<Mat2X> v, Ref<Mat2X> dq, Ref<Mat2X> dv, BlockTriMat* Kq, BlockTriMat* Kv) const {

   // Compute forces
   assemble(q, v, dv, Kq, Kv); 

   // Set change in positions based on velocity
   dq = v; 

   // Set change in velocities based on forces
   // TODO; we should probably precompute minv! 
   Eigen::Map<const Eigen::RowVectorXd> m(m_.data(), m_.size());
   Vec minv = m.cwiseInverse();

   dv.array().rowwise() /= m.array();

   if (Kq) Kq->multiply_rowwise( minv );
   if (Kv) Kv->multiply_rowwise( minv );

}

double LineModel::gravitational_energy(CRef<Mat2X> q) const
{
   Eigen::Map<const Eigen::VectorXd> m(m_.data(), m_.size());
   return q.row(1).dot(m)*params_.gravity.norm();
}
double LineModel::elastic_energy(CRef<Mat2X> q) const
{
   const Eigen::Index num_edges = line_.num_elements();

   double out = 0.0;

   for (Eigen::Index e = 0; e < num_edges; ++e){
      const Vec2 dq = edge_vector(q, e);
      const double len = dq.norm();

      const double s_main = len - line_.l_main[e];
      const double s_backup = len - line_.l_backup[e];

      if (s_main > 0) out += s_main*s_main*line_.k_main[e];
      if (s_backup > 0) out += s_backup*s_backup*line_.k_backup[e];
   }

   return 0.5*out;

   
}

double LineModel::kinetic_energy(CRef<Mat2X> v) const
{
   Eigen::Map<const Eigen::VectorXd> m(m_.data(), m_.size());
   
   return 0.5 * m.dot(v.colwise().squaredNorm());
}
