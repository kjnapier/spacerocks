#include "spacerocks.hpp"

struct Vector3 calc_vovec_from_kep(double mu, double a, double e, double r, double E){

  struct Vector3 vec;
  if (e < 1) {
    double xi = sqrt(mu * a) / r;
    vec.x = -xi * sin(E);
    vec.y = xi * sqrt(1 - e*e) * cos(E);
    vec.z = 0;
  }

  else {
    double xi = sqrt(-mu * a) / r;
    vec.x = -xi * sinh(E);
    vec.y = xi * sqrt(e*e - 1) * cosh(E);
    vec.z = 0;
  }

  return vec;

}