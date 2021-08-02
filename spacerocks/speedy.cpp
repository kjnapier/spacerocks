
const double mu_bary = 0.00029630927493457475;

struct StateVector{
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
};



//double calc_E(double e, double M){
//    /// Compute Danby (1988) fourth-order root-finding method with given number of steps
//    // This has quartic convergence.
//    // The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)
//
//    double k = 0.85;
//    double f_E, fP_E, fPP_E, fPPP_E, this_ell, old_E, delta_i1, delta_i2, delta_i3, esinE, ecosE;
//
//    // Define initial estimate
//    if((sin(M)) < 0) old_E = this_ell - k*e;
//    else old_E = this_ell + k*e;
//
//    // Perform Newton-Raphson estimate
//    for(int j = 0; j < 3; j++) {
//
//      // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
//      esinE = e*sin(old_E);
//      ecosE = e*cos(old_E);
//      f_E = old_E - esinE-this_ell;
//      fP_E = 1. - ecosE;
//      fPP_E = esinE;
//      fPPP_E = ecosE;
//
//      delta_i1 = -f_E/fP_E;
//      delta_i2 = -f_E/(fP_E+1./2.*delta_i1*fPP_E);
//      delta_i3 = -f_E/(fP_E+1./2.*delta_i2*fPP_E+1./6.*fPPP_E*delta_i2*delta_i2);
//
//      // Update E
//      old_E += delta_i3;
//    }
//
//    // Add to array
//    return old_E;
//}

double calc_E(double e, double M) {

  return E
}

StateVector kep_to_xyz_temp_cpp(double a, double e, double inc, double arg, double node, double M) {

  double E, true_anomaly, r, c, ox, oy, vox, voy
  double si, sa, sn, ci, ca, cn
  double c1, c2, c3, c4, c5, c6
  double x, y, z, vx, vy, vz


  if (e < 1) {

    // Compute E

    true_anomaly = 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    r = a * (1 - e * cos(E));

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(true_anomaly);
    oy = r * sin(true_anomaly);
    ovx = - c * sin(E);
    ovy = c * sqrt(abs(1 - e**2)) * cos(E);

  } else if (e >= 1) {

    // Compute E

    true_anomaly = 2 * arctan2(sqrt(e + 1) * tanh(E / 2), sqrt(e - 1));
    r = abs(a) * abs(1 - e * e) / (1 + e * cos(true_anomaly));

    c = sqrt(mu_bary * abs(a)) / r;

    ox = r * cos(true_anomaly);
    oy = r * sin(true_anomaly)
    ovx = - c * sinh(E);
    ovy = c * sqrt(abs(1 - e * e)) * cosh(E);

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

  x = ox * c1 - oy * c2;
  y = ox * c3 + oy * c4;
  z = ox * c5 + oy * c6;

  vx = vox * c1 - voy * c2;
  vy = vox * c3 + voy * c4;
  vz = vox * c5 + voy * c6;

  StateVector result; // = {x, y, z, vx, vy, vz};
  result.x = x;
  result.y = y;
  result.z = z;
  result.vx = vx;
  result.vy = vy;
  result.vz = vz;

  return result;

}

extern "C" {
    StateVector kep_to_xyz_temp(double a, double e, double inc, double arg, double node, double M)
    {
        return kep_to_xyz_temp_cpp(a, e, inc, arg, node, M);
    }
}
