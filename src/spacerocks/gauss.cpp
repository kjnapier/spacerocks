//#include "spacerocks.hpp"

#include <Eigen/Dense>
//#include <unsupported/Eigen/Polynomials>

void gauss() {

    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(0, 0, 0);
    Eigen::Matrix<double, 3, 3> c;
    // R2, R3;
    // R1 << 0.0, 0.0, 0.0;
    // R2 << 0.0, 0.0, 0.0;
    // R3 << 0.0, 0.0, 0.0;
    // Eigen::Vector3d R2(0, 0, 0); //detections[1].observer_position;
    // Eigen::Vector3d R3(0, 0, 0); //detections[2].observer_position;

    // Eigen::Vector3d rho1(0, 0, 0); //detections[0].pointing_vector;
    // Eigen::Vector3d rho2(0, 0, 0); //detections[1].pointing_vector;
    // Eigen::Vector3d rho3(0, 0, 0); //detections[2].pointing_vector;

    // Eigen::Vector3d p1 = rho2.cross(rho3);
    // Eigen::Vector3d p2 = rho1.cross(rho3);
    // Eigen::Vector3d p3 = rho1.cross(rho2);

    // double D0 = rho1.dot(p1);

    // Eigen::Matrix<double, 3, 3> D;
    // D(0, 0) = R1.dot(p1);
    // D(0, 1) = R1.dot(p2);
    // D(0, 2) = R1.dot(p3); 
    // D(1, 0) = R2.dot(p1);
    // D(1, 1) = R2.dot(p2);
    // D(1, 2) = R2.dot(p3);
    // D(2, 0) = R3.dot(p1);
    // D(2, 1) = R3.dot(p2);
    // D(2, 2) = R3.dot(p3);

    // double A = (1/D0) * (-D(0, 1) * (tau3/tau) + D(1, 1) + D(2, 1) * (tau1/tau));
    // double B = (1/(6 * D0)) * (D(0, 1) * (tau3**2 - tau**2) * (tau3/tau) + D(2, 1) * (tau**2 - tau1**2) * (tau1/tau));
    // double E = R2.dot(rho2);

    // double R2sq = R2.dot(R2);

    // double a = -(A * A + 2 * A * E + R2sq);
    // double b = -2 * mu_bary * B * (A + E);
    // double c = -mu_bary * mu_bary * B * B;

    // Eigen::Matrix<double, 8, 1> coeff;
    // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    // coeff << c, 0.0, 0.0, b, 0.0, a, 0.0, 1.0;
    // solver.compute(coeff);
    // const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();

    // for (auto& root : r) {
    //     if (root.imag() == 0.0) {
    //         double a1 = (1/D0) * ((6 * (D(2, 0) * (tau1/tau3) + D(1, 0) * (tau/tau3)) * root.real()**3 + mu_bary * D(2, 0) * (tau**2 - tau1**2) * (tau1/tau3)) / (6 * root.real()**3 + mu_bary * (tau**2 - tau3**2)) - D(0, 0));
    //         double a2 = A + (mu_bary * B) / root.real()**3;
    //         double a3 = (1/D0) * ((6 * (D(0, 2) * (tau3/tau1) - D(1, 2) * (tau/tau1)) * root.real()**3 + mu_bary * D(0, 2) * (tau**2 - tau3**2) * (tau3/tau1)) / (6 * root.real()**3 + mu_bary * (tau**2 - tau1**2)) - D(2, 2));

    //         double f1 = 1 - 0.5 * (mu_bary/root.real()**3) * tau1**2;
    //         double f3 = 1 - 0.5 * (mu_bary/root.real()**3) * tau3**2;
    //         double g1 = tau1 - (1/6) * (mu_bary / root.real()**3) * tau1**3;
    //         double g3 = tau3 - (1/6) * (mu_bary / root.real()**3) * tau3**3;

    //         Eigen::Vector3d r1 = R1 + a1 * rho1;
    //         Eigen::Vector3d r2 = R2 + a2 * rho2;
    //         Eigen::Vector3d r3 = R3 + a3 * rho3;

    //         Eigen::Vector3d v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1);

    //         x, y, z = r2(0), r2(1), r2(2);
    //         vx, vy, vz = v2(0), v2(1), v2(2);

    //         StateVector result = {x, y, z, vx, vy, vz};
    //         results.push_back(result);
    //     }
    // }
    // return results;
    return;
}

int main() {
    gauss();
    return 0;
}