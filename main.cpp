#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "quaternion.hpp"
#include "utils.hpp"

using namespace Eigen;

constexpr double MU = 4902.800238e9;   // m^3/s^2
constexpr double RMOON = 1727.4e3;     // m
constexpr double dt = 0.05;            // sec
constexpr double Tfinal = 100.0;       // sec
constexpr int N = static_cast<int>(Tfinal / dt);

Vector3d WMOON(0, 0, 2.597216148801086e-6);
Matrix3d I3 = Matrix3d::Identity();

void computeFandPhi(const Vector3d& r, const Vector3d& v, const Vector4d& q, const Vector3d& wb, const Vector3d& ab,
                    MatrixXd& F, MatrixXd& Phi) 
{
    using namespace Eigen;

    const double MU = 4902.800238e9; // m^3/s^2
    Matrix3d I3 = Matrix3d::Identity();
    F = MatrixXd::Zero(24, 24);

    // Partial gravity gradient
    Matrix3d g_rr = (3 * MU / pow(r.norm(), 5)) * (r * r.transpose()) - (MU / pow(r.norm(), 3)) * I3;

    // Skew symmetric matrices
    Matrix3d Omega_wb = cross_prod(wb);
    Matrix3d Omega_ab = cross_prod(ab);

    // Partial derivatives of inertial accel wrt theta (attitude error)
    Matrix3d A_th = qtom(qconj(q)) * cross_prod(ab);

    // Populate F
    F.block<3,3>(0,3) = I3;        // r_dot = v
    F.block<3,3>(3,0) = g_rr;      // dv/dr
    F.block<3,3>(3,6) = A_th;      // dv/dtheta
    F.block<3,3>(3,21) = qtom(qconj(q)); // dv/dba
    F.block<3,3>(6,6) = -cross_prod(wb); // dtheta/dtheta
    F.block<3,3>(6,10) = -I3;      // dtheta/dbg
    F.block<3,3>(6,9) = -Omega_wb; // dtheta/dsg
    F.block<3,3>(3,15) = qtom(qconj(q)) * Omega_ab; // dv/dsa
    F.block<3,3>(3,12) = -qtom(qconj(q)) * cross_prod(ab); // dv/deo

    // Discretization using matrix exponential
    // Phi = (F * dt).exp();  // Requires Eigen unsupported module: MatrixFunctions

    // Alternative if exp() is unavailable:
    Phi = MatrixXd::Identity(24,24) + F * dt + 0.5 * F * F * dt * dt;
}

void runCodeErrorModels(MatrixXd& imuData)
{
    std::cout << "Running code error models using linear covariance analysis" << std::endl;
    std::cout << "Output: (Position Error NED (m), Velocity Error NED (mps))" << std::endl;

    std::vector<std::string> labels = {
        "gyro scale factor", "gyro bias", "gyro random walk",
        "accel misalignment", "accel scale factor", "accel bias",
        "accel random walk", "gyro rss", "accel rss"
    };

    MatrixXd toggle = MatrixXd::Zero(7, 9);

    toggle.block<7,7>(0,0) = MatrixXd::Identity(7,7);
    toggle.col(7) << 1, 1, 1, 0, 0, 0, 0;
    toggle.col(8) << 1, 1, 1, 1, 1, 1, 1;

    for (int kk = 0; kk < 9; ++kk) {
        // Initial position/velocity
        Vector3d r(RMOON, 0, 0);
        Vector3d v = WMOON.cross(r);

        // Initial attitude quaternion
        Vector4d q(0, 0, 0, 1);

        std::vector<double> tsav(N);
        MatrixXd sigma(6, N);

        // Initial P
        VectorXd r_var = VectorXd::Zero(3);
        VectorXd v_var = VectorXd::Zero(3);
        VectorXd theta_var = VectorXd::Zero(3);

        VectorXd sg_var = toggle(0,kk) * pow(5.0 / 1e6, 2) * VectorXd::Ones(3);
        VectorXd bg_var = toggle(1,kk) * pow(deg2rad(0.05 / 60.0 / 60.0), 2) * VectorXd::Ones(3);
        double v_gyro = toggle(2,kk) * pow(deg2rad(0.0001), 2);

        VectorXd eo_var = toggle(3,kk) * pow(deg2rad(15.0 / 60.0 / 60.0), 2) * VectorXd::Ones(3);
        VectorXd sa_var = toggle(4,kk) * pow(175.0 / 1e6, 2) * VectorXd::Ones(3);
        VectorXd ba_var = toggle(5,kk) * pow(100.0 * 9.81e-6, 2) * VectorXd::Ones(3);
        double v_accel = toggle(6,kk) * pow(0.0005 * 0.3048, 2);

        MatrixXd G = MatrixXd::Zero(24, 6);
        G.block<3,3>(3, 0) = I3;
        G.block<3,3>(6, 3) = I3;

        MatrixXd Q = MatrixXd::Zero(6, 6);
        Q.block<3,3>(0, 0) = v_accel * I3;
        Q.block<3,3>(3, 3) = v_gyro * I3;

        VectorXd diag_P(24);
        diag_P << r_var, v_var, theta_var, sg_var, bg_var, eo_var, sa_var, ba_var;
        MatrixXd P = diag_P.asDiagonal();

        double t = 0.0;
        for (int i = 0; i < N - 1; ++i) {
            Vector3d wb = Vector3d::Zero();
            Vector3d ab1 = Vector3d::Zero();
            Vector3d ab2 = Vector3d::Zero();

            if (t <= 30.0) {
                wb = imuData.row(i).segment<3>(1).transpose();
                ab1 = imuData.row(i).segment<3>(4).transpose();
                ab2 = imuData.row(i+1).segment<3>(4).transpose();
                // Add errors here if needed
            }

            tsav[i] = t;

            // Attitude propagation
            Vector4d qdt = qnorm(q + qdot(q, wb) * dt);

            // Gravity
            Vector3d g = -MU * r / pow(r.norm(), 3);

            // Acceleration
            Vector3d a_avg = 0.5 * qvq(qconj(qdt), ab2) + 0.5 * qvq(qconj(q), ab1) + g;

            // Position and velocity propagation
            r += v * dt + 0.5 * a_avg * dt * dt;
            v += a_avg * dt;
            q = qdt;

            // Propagate covariance matrix
            MatrixXd F, Phi;
            computeFandPhi(r, v, q, wb, ab1, F, Phi);
            P = Phi * P * Phi.transpose() + G * Q * G.transpose() * dt;

            Matrix3d Tf = inertialToNED(r);
            MatrixXd A = MatrixXd::Zero(3, 24);
            A.block<3,3>(0,0) = Tf;
            MatrixXd P_r_ned = A * P * A.transpose();
            MatrixXd B = MatrixXd::Zero(3, 24);
            B.block<3,3>(0,3) = Tf;
            MatrixXd P_v_ned = B * P * B.transpose();

            sigma.col(i) << P_r_ned.diagonal().cwiseSqrt(), P_v_ned.diagonal().cwiseSqrt();
            t += dt;
        }

        std::cout << "Final Sigma for " << labels[kk] << ":   " << sigma.col(N-2).transpose() << "\n";

        //     // Print final state
        // std::cout.precision(12);
        // std::cout << "Final position:\n" << r.transpose() << std::endl;
        // std::cout << "Final velocity:\n" << v.transpose() << std::endl;
        // std::cout << "Final quaternion:\n" << q.transpose() << std::endl;

        // // Expected final values (no errors)
        // Vector3d rf(1.85560674378386e6, 0.14361820808053e6, 0);
        // Vector3d vf(1.40261326327722e3, 1.81913089494382e3, 0);
        // Vector4d qf(0, 0, 0.70710825356386, 0.70710530880617);

        // std::cout << "\nPosition Error: " << (r - rf).transpose() << std::endl;
        // std::cout << "Velocity Error: " << (v - vf).transpose() << std::endl;
        // std::cout << "Quaternion Error: " << (q - qf).transpose() << std::endl;
    }
}

void runSimulation(MatrixXd& imuData)
{
    std::cout << "Running Rocket Simulation" << std::endl;
    // Initialize state
    Vector3d r(RMOON, 0, 0); // Position
    Vector3d v = WMOON.cross(r); // Velocity

    Vector4d q(0, 0, 0, 1); // Quaternion
    q = qnorm(q);

    double t = 0;

    for (int i = 0; i < N - 1; ++i) {
        Vector3d wb = Vector3d::Zero();
        Vector3d ab1 = Vector3d::Zero();
        Vector3d ab2 = Vector3d::Zero();


        if (t <= 30.0 && i + 1 < imuData.size()) {
           wb = imuData.row(i).segment<3>(1).transpose();
           ab1 = imuData.row(i).segment<3>(4).transpose();
           ab2 = imuData.row(i+1).segment<3>(4).transpose();
        }

        // Propagate attitude
        Vector4d qdt = qnorm(q + qdot(q, wb) * dt);

        // Gravity
        Vector3d g = -MU * r / pow(r.norm(), 3);

        // Acceleration
        Vector3d a_avg = 0.5 * qvq(qconj(qdt), ab2) + 0.5 * qvq(qconj(q), ab1) + g;

        // Update position and velocity
        r += v * dt + 0.5 * a_avg * dt * dt;
        v += a_avg * dt;
        q = qdt;

        t += dt;
    }

    // Print final state
    std::cout.precision(12);
    std::cout << "Final position:\n" << r.transpose() << std::endl;
    std::cout << "Final velocity:\n" << v.transpose() << std::endl;
    std::cout << "Final quaternion:\n" << q.transpose() << std::endl;

    // Expected final values (no errors)
    Vector3d rf(1.85560674378386e6, 0.14361820808053e6, 0);
    Vector3d vf(1.40261326327722e3, 1.81913089494382e3, 0);
    Vector4d qf(0, 0, 0.70710825356386, 0.70710530880617);

    std::cout << "\nPosition Error: " << (r - rf).transpose() << std::endl;
    std::cout << "Velocity Error: " << (v - vf).transpose() << std::endl;
    std::cout << "Quaternion Error: " << (q - qf).transpose() << std::endl;
}

int main() 
{
    MatrixXd imuData = loadIMUData("IMUdata.csv");

    runCodeErrorModels(imuData);

    std::cout << std::endl << std::endl;

    runSimulation(imuData);

    return 0;
}
