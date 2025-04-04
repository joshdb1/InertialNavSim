#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace Eigen;
using namespace std;

// Constants
const double MU = 4902.800238e9;
const double RMOON = 1727.4e3;
const Vector3d WMOON(0, 0, 2.597216148801086e-6);
const double dt = 0.05;
const double Tfinal = 100;
const int N = Tfinal / dt;

// Quaternion utilities
Vector4d qnorm(const Vector4d &q) {
    return q / q.norm();
}

Vector4d qdot(const Vector4d &q, const Vector3d &w) {
    Matrix4d Omega;
    Omega <<   0,   -w(0), -w(1), -w(2),
               w(0),  0,    w(2), -w(1),
               w(1), -w(2),  0,    w(0),
               w(2),  w(1), -w(0),  0;
    return 0.5 * Omega * q;
}

Vector4d qconj(const Vector4d &q) {
    return Vector4d(-q(0), -q(1), -q(2), q(3));
}

Vector3d qvq(const Vector4d &q, const Vector3d &v) {
    Quaterniond Q(q(3), q(0), q(1), q(2)); // Eigen uses w,x,y,z
    return Q * v;
}

// Load CSV data
vector<VectorXd> loadIMUdata(const string &filename) {
    vector<VectorXd> data;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        istringstream ss(line);
        string val;
        VectorXd row(7);
        for (int i = 0; i < 7; ++i) {
            getline(ss, val, ',');
            row(i) = stod(val);
        }
        data.push_back(row);
    }
    return data;
}

int main() {
    // Load IMU data
    vector<VectorXd> IMUdata = loadIMUdata("IMUdata.csv");

    // Initialize state
    Vector3d r(RMOON, 0, 0); // Position
    Vector3d v = WMOON.cross(r); // Velocity

    Vector4d q(0, 0, 0, 1); // Quaternion
    q = qnorm(q);

    // Storage for plots
    vector<double> tsav, hsav, vsav, gamsav, gsav;

    double t = 0;

    for (int i = 0; i < N - 1; ++i) {
        Vector3d wb = Vector3d::Zero(), ab1 = Vector3d::Zero(), ab2 = Vector3d::Zero();

        if (t <= 30.0 && i + 1 < IMUdata.size()) {
            wb = IMUdata[i].segment<3>(1);
            ab1 = IMUdata[i].segment<3>(4);
            ab2 = IMUdata[i + 1].segment<3>(4);

            // Insert gyro/accel error models here (currently zero)
        }

        // Save variables
        tsav.push_back(t);
        hsav.push_back(r.norm() - RMOON);
        vsav.push_back(v.norm());
        gamsav.push_back(atan2(v(0), v(1)));
        gsav.push_back(ab1.norm() / 9.8);

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

    // Plotting
    plt::subplot(2, 2, 1);
    plt::plot(tsav, hsav);
    plt::xlabel("Time (sec)");
    plt::ylabel("Altitude (m)");

    plt::subplot(2, 2, 2);
    plt::plot(tsav, vsav);
    plt::xlabel("Time (sec)");
    plt::ylabel("Speed (m)");

    plt::subplot(2, 2, 3);
    plt::plot(tsav, gsav);
    plt::xlabel("Time (sec)");
    plt::ylabel("G-Load (g)");

    plt::subplot(2, 2, 4);
    vector<double> gam_deg(gamsav.size());
    transform(gamsav.begin(), gamsav.end(), gam_deg.begin(), [](double rad){ return rad * 180.0 / M_PI; });
    plt::plot(tsav, gam_deg);
    plt::xlabel("Time (sec)");
    plt::ylabel("Flight Path Angle (deg)");

    plt::show();

    // Print final state
    cout.precision(12);
    cout << "Final position:\n" << r.transpose() << endl;
    cout << "Final velocity:\n" << v.transpose() << endl;
    cout << "Final quaternion:\n" << q.transpose() << endl;

    // Expected final values (no errors)
    Vector3d rf(1.85560674378386e6, 0.14361820808053e6, 0);
    Vector3d vf(1.40261326327722e3, 1.81913089494382e3, 0);
    Vector4d qf(0, 0, 0.70710825356386, 0.70710530880617);

    cout << "\nPosition Error: " << (r - rf).transpose() << endl;
    cout << "Velocity Error: " << (v - vf).transpose() << endl;
    cout << "Quaternion Error: " << (q - qf).transpose() << endl;

    return 0;
}
