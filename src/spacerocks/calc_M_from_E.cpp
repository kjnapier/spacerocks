#include "spacerocks.hpp"

double calc_M_from_E(double e, double E) {

  double M;


  if (e < 1) {
    M = E - e * sin(E);
  }
  else {
    M = e * sinh(E) - E;
  }

  return M;

}
