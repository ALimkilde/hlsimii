#pragma once

#include <cassert>
#include <Eigen/Dense>


using Vec = Eigen::VectorXd;
using Vec2 = Eigen::Vector2d;

using Mat = Eigen::MatrixXd;
using Mat2X = Eigen::Matrix2Xd;
using Mat2 = Eigen::Matrix2d;


template <typename T> using Ref  = Eigen::Ref<T>;
template <typename T> using CRef = Eigen::Ref<const T>;


// Views of a flat vector [x0 y0 x1 y1 ...] as a 2 x n matrix (no copy)
inline Eigen::Map<Mat2X> as_mat2x(Ref<Vec> x) {
   assert(x.size() % 2 == 0);
   return Eigen::Map<Mat2X>(x.data(), 2, x.size() / 2);
}

inline Eigen::Map<const Mat2X> as_const_mat2x(CRef<Vec> x) {
   assert(x.size() % 2 == 0);
   return Eigen::Map<const Mat2X>(x.data(), 2, x.size() / 2);
}
