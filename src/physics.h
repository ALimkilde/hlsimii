#pragma once

#include <Eigen/Dense>

#include <line.h>
#include <mesh.h>

using Vec2 = Eigen::Vector2d;
using Vec = Eigen::VectorXd;

struct Params {
   Vec2 gravity{0.0, -9.81};
   // Drag and dampening parameters
};

struct workspace {
   Vec state; // [x1 y1 x2 y2 ...| vx1 vy1 vx2 vy2 ...]
};

class LineModel {
  public:
      explicit LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p);  // sizes all buffers here, once
  
      Eigen::Index size() const;
 
      void rhs(double t, const Vec& y, Vec& dy); // y input state; dy output rhs;
  private: 
      DiscreteLine line_;
      Vec2 anchor_left_, anchor_right_;
      Params params_;                                               
      Eigen::Matrix2Xd Fe_;  // per-edge Forces
};
