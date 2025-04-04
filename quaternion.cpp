#include "quaternion.hpp"
#include "utils.hpp"

// Quaternion multiplication
Vector4d qmult(const Vector4d &q1, const Vector4d &q2) 
{
    Vector4d result;
    result(0) = q1(3)*q2(0) + q1(0)*q2(3) + q1(1)*q2(2) - q1(2)*q2(1);
    result(1) = q1(3)*q2(1) - q1(0)*q2(2) + q1(1)*q2(3) + q1(2)*q2(0);
    result(2) = q1(3)*q2(2) + q1(0)*q2(1) - q1(1)*q2(0) + q1(2)*q2(3);
    result(3) = q1(3)*q2(3) - q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2);
    return result;
}

Vector4d qxq(const Vector4d& a, const Vector4d& b) 
{
    // Create the 4x4 matrix B based on quaternion b
    Matrix4d B;
    B << b(3),  b(2), -b(1),  b(0),
        -b(2),  b(3),  b(0),  b(1),
         b(1), -b(0),  b(3),  b(2),
        -b(0), -b(1), -b(2),  b(3);
    
    // Perform the matrix-vector multiplication
    Vector4d c = B * a;
    
    return c;
}

// Quaternion derivative
Vector4d qdot(const Vector4d& q, const Vector3d& omega) 
{
    Vector4d qw;
    qw << omega, 0;
    return 0.5 * qxq(q, qw);
}

// Quaternion normalization
Vector4d qnorm(const Vector4d &q) 
{
    return q / q.norm();
}

// Quaternion conjugate
Vector4d qconj(const Vector4d &q) 
{
    Vector4d qc;
    qc = -q;
    qc(3) = -q(3);
    return qc;
}

// Rotate vector v by quaternion q
Vector3d qvq(const Vector4d& q, const Vector3d& v) 
{
    // Convert v to quaternion with 0 scalar part
    Vector4d v_quat;
    v_quat << v, 0;

    // Conjugate of q
    Vector4d q_conj;
    q_conj << -q.head<3>(), q(3);

    // Quaternion multiplication: q * v_quat
    Vector4d temp;
    Vector3d qv = q.head<3>();
    double qw = q(3);
    Vector3d temp_v = qw * v + qv.cross(v);
    double temp_w = -qv.dot(v);
    temp << temp_v, temp_w;

    // Now temp * q_conj
    Vector3d qcv = q_conj.head<3>();
    double qcw = q_conj(3);
    Vector3d out_v = temp_w * qcv + qcw * temp.head<3>() + temp.head<3>().cross(qcv);

    return out_v;
}

// Convert quaternion to rotation matrix
Matrix3d qtom(const Vector4d& q) 
{
    // Extract the vector part (qv) and scalar part (q4)
    Eigen::Vector3d qv = q.head<3>();  // The vector part of the quaternion
    double q4 = q(3);  // The scalar part of the quaternion
    
    // Calculate cross product matrix of qv
    Eigen::Matrix3d Q = cross_prod(qv);

    // Calculate the matrix m
    double qv_norm_sq = qv.squaredNorm();  // qv' * qv (norm squared)
    double q4_squared = q4 * q4;

    Eigen::Matrix3d m = (q4_squared - qv_norm_sq) * Eigen::Matrix3d::Identity()  // Diagonal term
                        + 2 * qv * qv.transpose()  // Outer product (2 * qv * qv')
                        - 2 * q4 * Q;  // Cross product matrix term (-2 * q(4) * Q)
    
    return m;
}