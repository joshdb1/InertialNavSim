#include "utils.hpp"
#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>

Eigen::MatrixXd loadIMUData(const std::string& filePath) 
{
    // Open the file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening IMU data file!" << std::endl;
        exit(EXIT_FAILURE); // Exit if the file can't be opened
    }

    std::vector<std::vector<double>> data;
    std::string line;
    
    // Read each line of the file
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        std::string value;
        
        // Read the values in each line (assuming comma-separated values)
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }
        
        data.push_back(row);
    }
    
    // Close the file
    file.close();
    
    // Convert to Eigen matrix
    Eigen::MatrixXd imuData(data.size(), data[0].size());
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            imuData(i, j) = data[i][j];
        }
    }

    return imuData; // Return the loaded data as an Eigen matrix
}

Eigen::Matrix3d cross_prod(const Eigen::Vector3d& v) 
{
    Eigen::Matrix3d S;
    S <<     0, -v(2),  v(1),
          v(2),     0, -v(0),
         -v(1),  v(0),     0;
    return S;
}

// Function to compute the rotation matrix from Inertial to NED frame
Matrix3d inertialToNED(const Vector3d& r) 
{
    // Normalize the inertial position vector
    Vector3d id = -r.normalized();  // Inertial to Earth frame (downward direction)

    // Define the unit vector along the z-axis (upward)
    Vector3d iz(0, 0, 1);

    // Compute the East unit vector (cross product of id and iz)
    Vector3d E = id.cross(iz);
    Vector3d ie = E.normalized();  // Normalize the East vector

    // Compute the North unit vector (cross product of ie and id)
    Vector3d in = ie.cross(id);

    // Construct the transformation matrix from Inertial to NED frame
    Matrix3d T;
    T << in, ie, id;

    return T;
}

double deg2rad(double degrees) 
{
    return degrees * (M_PI / 180.0);
}