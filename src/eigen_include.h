#pragma once

#include <Eigen/Dense>


using Vec = Eigen::VectorXd;
using Vec2 = Eigen::Vector2d;

using Mat = Eigen::MatrixXd;
using Mat2X = Eigen::Matrix2Xd;
using Mat2 = Eigen::Matrix2d;


template <typename T> using Ref  = Eigen::Ref<T>;
template <typename T> using CRef = Eigen::Ref<const T>;
