#pragma once

#include <Eigen/Dense>

using namespace Eigen;

Eigen::MatrixXd loadIMUData(const std::string& filePath);

Eigen::Matrix3d cross_prod(const Eigen::Vector3d& v);

Matrix3d inertialToNED(const Vector3d& r);

double deg2rad(double degrees);