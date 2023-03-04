#include "spacerocks.hpp"

double calc_E_from_f(double e, double f) {

  double E;

  if (e < 1) {
    E = 2 * atan2(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2));
  }
  else {
    double cta = cos(f);
    E = acosh((cta + e) / (1 + e * cta));
    if (f < 0) {
      E *= -1;
    }
  }
  return E;
}



