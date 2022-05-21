#include "spacerocks.hpp"

double calc_f_from_E(double e, double E) {

  double f;

  if (e < 1) {
    f = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    //f = 2 * atan2(sqrt((1 + e)/(1 - e)) * tan(E / 2), 1);
  }
  else {
    //f = 2 * atan2(sqrt(e + 1) * sinh(E / 2), sqrt(e - 1) * cosh(E / 2));
    f = 2 * atan2(sqrt((e + 1)/(e - 1)) * tanh(E / 2), 1);
  }

  return f; 

}