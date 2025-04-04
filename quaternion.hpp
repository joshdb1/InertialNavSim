// quaternion.hpp
#pragma once

#include <Eigen/Dense>

using namespace Eigen;

// Quaternion multiplication
Vector4d qmult(const Vector4d &q1, const Vector4d &q2);

// Quaternion cross product
Vector4d qxq(const Vector4d& a, const Vector4d& b); 

// Quaternion derivative
Vector4d qdot(const Vector4d& q, const Vector3d& omega); 

// Quaternion normalization
Vector4d qnorm(const Vector4d &q); 

// Quaternion conjugate
Vector4d qconj(const Vector4d &q); 

// Rotate vector v by quaternion q
Vector3d qvq(const Vector4d& q, const Vector3d& v); 

// Convert quaternion to rotation matrix
Matrix3d qtom(const Vector4d& q); 