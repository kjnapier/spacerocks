#include "spacerocks.hpp"

double compute_lst(double epoch, double lon) {

    double T, theta;

    T = (epoch - 2451545.0) / 36525;
    theta = 280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0);
    return theta + lon;

}