#include <math.h>
#include <stdlib.h>

const double mu_bary = 0.00029630927493457475;

struct StateVector{
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
};

double calc_E(double e, double M){
    /// Compute Danby (1988) fourth-order root-finding method with given number of steps
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
      }
    }

    else {

      E = M/fabs(M)*log(2.*fabs(M)/e + 1.8);

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


struct StateVector kep_to_xyz_temp_cpp(double a, double e, double inc, double arg, double node, double M) {

  double E, true_anomaly, r, c, ox, oy, vox, voy;
  double cosE, aba, omece;
  double si, sa, sn, ci, ca, cn;
  double c1, c2, c3, c4, c5, c6;

  struct StateVector rock;

  if (e < 1) {

    // Compute E
    E = calc_E(e, M);


    cosE = cos(E);
    omece = 1 - e * cosE;

    //true_anomaly = acos((cosE - e) / omece);
    true_anomaly = 2 * atan2(sqrt((1 + e)/(1 - e)) * sin(E / 2), cos(E / 2));
    r = a * omece;

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(true_anomaly);
    oy = r * sin(true_anomaly);
    vox = - c * sin(E);
    voy = c * sqrt(1 - e * e) * cosE;

  } else {

    aba = fabs(a);

    E = calc_E(e, M);

    true_anomaly = 2 * atan2(sqrt(e + 1) * tanh(E / 2), sqrt(e - 1));
    r = aba * fabs(1 - e * e) / (1 + e * cos(true_anomaly));

    c = sqrt(mu_bary * aba) / r;

    ox = r * cos(true_anomaly);
    oy = r * sin(true_anomaly);
    vox = - c * sinh(E);
    voy = c * sqrt(fabs(1 - e * e)) * cosh(E);

  }

  sa = sin(arg);
  si = sin(inc);
  sn = sin(node);
  ca = cos(arg);
  ci = cos(inc);
  cn = cos(node);

  c1 = ca * cn - sa * sn * ci;
  c2 = sa * cn + ca * sn * ci;
  c3 = ca * sn + sa * cn * ci;
  c4 = ca * cn * ci - sa * sn;
  c5 = sa * si;
  c6 = ca * si;

  rock.x = ox * c1 - oy * c2;
  rock.y = ox * c3 + oy * c4;
  rock.z = ox * c5 + oy * c6;
  rock.vx = vox * c1 - voy * c2;
  rock.vy = vox * c3 + voy * c4;
  rock.vz = vox * c5 + voy * c6;

  return rock;

}


double* kep_to_xyz_temp(int N, double *as, double *es, double *incs, double *args, double *nodes, double *Ms)
{

  double* output = malloc(N * 6 * sizeof(double));
  int dummy;

  struct StateVector rock;

  for (int idx = 0; idx < N; idx++) {

    rock = kep_to_xyz_temp_cpp(as[idx], es[idx], incs[idx], args[idx], nodes[idx], Ms[idx]);

    dummy = idx * 6;

    output[dummy] = rock.x;
    output[dummy + 1] = rock.y;
    output[dummy + 2] = rock.z;
    output[dummy + 3] = rock.vx;
    output[dummy + 4] = rock.vy;
    output[dummy + 5] = rock.vz;

  }

  return output;

}
