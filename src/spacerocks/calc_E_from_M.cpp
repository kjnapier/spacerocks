#include "spacerocks.hpp"

double calc_E_from_M(double e, double M){
    // Compute Danby (1988) fourth-order root-finding method with given number of steps
    // This has quartic convergence.
    // The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)

    double E;

    if (e < 1) {

      double k = 0.85;
      double f_E, fP_E, fPP_E, fPPP_E, delta_i1, delta_i2, delta_i3, esinE, ecosE;

      // Define initial estimate
      if((sin(M)) < 0) E = M - k * e;
      else E = M + k * e;

      // Perform Newton-Raphson estimate
      for(int j = 0; j < 10; j++) {

        // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
        esinE = e*sin(E);
        ecosE = e*cos(E);
        f_E = E - esinE - M;
        fP_E = 1 - ecosE;
        fPP_E = esinE;
        fPPP_E = ecosE;

        delta_i1 = -f_E/fP_E;
        delta_i2 = -f_E/(fP_E+1./2.*delta_i1*fPP_E);
        delta_i3 = -f_E/(fP_E+1./2.*delta_i2*fPP_E+1./6.*fPPP_E*delta_i2*delta_i2);
      
        // Update E
        E += delta_i3;
        if(fabs(delta_i3) < 1.e-16){
          break;
        }
      }
    }

    else {

      if (M != 0) {
        E = M/fabs(M)*log(2.*fabs(M)/e + 1.8);
      }
      else {
        E = 0;
      }
      
  		double F = E - e*sinh(E) + M;
  		for(int i=0; i < 100; i++){
  			E = E - F/(1.0 - e*cosh(E));
  			F = E - e*sinh(E) + M;
  			if(fabs(F) < 1.e-16){
  				break;
        }
      }
    }
    return E;
}