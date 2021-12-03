#include <math.h>
#include <stdlib.h>
#include <stdio.h>

const double EMIN = 1e-8;
const double IMIN = 1e-8;
const double mu_bary = 0.00029630927493457475;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

struct StateVector{
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
};

struct KeplerOrbit{
  double a;
  double e;
  double inc;
  double arg;
  double node;
  double f;
};

struct Vector3{
  double x;
  double y;
  double z;
};

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


struct StateVector kepM_to_xyz(double a, double e, double inc, double arg, double node, double M) {

  double E, f, r, c, ox, oy, vox, voy;
  double cosE, omece;
  double si, sa, sn, ci, ca, cn;
  double c1, c2, c3, c4, c5, c6;

  struct StateVector rock;

  if (e < 1) {

    // Compute E
    E = calc_E_from_M(e, M);


    cosE = cos(E);
    omece = 1 - e * cosE;

    //f = acos((cosE - e) / omece);
    f = 2 * atan2(sqrt((1 + e)/(1 - e)) * sin(E / 2), cos(E / 2));
    r = a * omece;

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sin(E);
    voy = c * sqrt(1 - e * e) * cosE;

  } else {

    E = calc_E_from_M(e, M);

    //f = 2 * atan2(sqrt(e + 1) * tanh(E / 2), sqrt(e - 1));
    f = 2 * atan2(sqrt(e + 1) * sinh(E / 2), cosh(E / 2) * sqrt(e - 1));
    r = a * (1 - e*e) / (1 + e * cos(f));

    c = sqrt(- mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sinh(E);
    voy = c * sqrt(e*e - 1) * cosh(E);

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

struct StateVector kepE_to_xyz(double a, double e, double inc, double arg, double node, double E) {

  double f, r, c, ox, oy, vox, voy;
  double cosE, omece;
  double si, sa, sn, ci, ca, cn;
  double c1, c2, c3, c4, c5, c6;

  struct StateVector rock;

  if (e < 1) {

    cosE = cos(E);
    omece = 1 - e * cosE;

    //f = acos((cosE - e) / omece);
    f = 2 * atan2(sqrt((1 + e)/(1 - e)) * sin(E / 2), cos(E / 2));
    r = a * omece;

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sin(E);
    voy = c * sqrt(1 - e * e) * cosE;

  } else {

    f = 2 * atan2(sqrt(e + 1) * tanh(E / 2), sqrt(e - 1));
    r = a * (1 - e*e) / (1 + e * cos(f));

    c = sqrt(- mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sinh(E);
    voy = c * sqrt(e*e - 1) * cosh(E);

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

struct KeplerOrbit calc_kep_from_xyz(double mu, double x, double y, double z, double vx, double vy, double vz) {

  struct KeplerOrbit kep;
  struct Vector3 position;
  struct Vector3 velocity;
  struct Vector3 evec;
  struct Vector3 nvec;
  struct Vector3 hvec;

  double r = sqrt(x*x + y*y + z*z);
  double vsq = vx*vx + vy*vy + vz*vz;
  double a, e, inc, arg, node, f;

  position.x = x;
  position.y = y;
  position.z = z;

  velocity.x = vx;
  velocity.y = vy;
  velocity.z = vz;

  // Compute the anguler momentum vector
  hvec.x = position.y * velocity.z - position.z * velocity.y;
  hvec.y = position.z * velocity.x - position.x * velocity.z;
  hvec.z = position.x * velocity.y - position.y * velocity.x;

  // Compute the eccentricity vector
  evec.x = (velocity.y * hvec.z - velocity.z * hvec.y) / mu - position.x / r;
  evec.y = (velocity.z * hvec.x - velocity.x * hvec.z) / mu - position.y / r;
  evec.z = (velocity.x * hvec.y - velocity.y * hvec.x) / mu - position.z / r;


  // Compute the normal vector
  nvec.x = -hvec.y;
  nvec.y = hvec.x;
  double n = sqrt(nvec.x*nvec.x + nvec.y*nvec.y);

  // Compute the semi-major axis (a)
  a = 1 / (2 / r - vsq / mu);

  // Compute the eccentricity (e)
  e = sqrt(evec.x*evec.x + evec.y*evec.y + evec.z*evec.z);

  // Compute the inclination (inc)
  inc = acos(hvec.z / sqrt(hvec.x*hvec.x + hvec.y*hvec.y + hvec.z*hvec.z));

  // Compute the longitude of ascending node (node)
  if (inc == 0) {
    node = 0;
  }
  else {
    node = acos(nvec.x / n);
  }
  if (nvec.y < 0) {
    node = 2 * M_PI - node;
  }

  // Compute the argument of pericenter (arg)
  if (e < EMIN) {
    arg = 0;
  }
  else if (inc < IMIN || inc > M_PI - IMIN) {
    double theta = acos(position.x / r);
    if (position.y < 0) {
      theta = 2 * M_PI - theta;
    }
    double varpi = acos(evec.x / e);
    if (evec.y < 0) {
      varpi = 2 * M_PI - varpi;
    }
    if (inc < M_PI / 2) {
      arg = varpi - node;
    }
    else {
      arg = node - varpi;
    }
  }
  else {
    arg = acos((nvec.x*evec.x + nvec.y*evec.y) / (n * e));
    if (evec.z < 0) {
      arg = 2 * M_PI - arg;
    }
  }

  // Compute the true anomaly (f)
  if (e < 1) {
    if (inc < IMIN || inc > M_PI - IMIN) {
      // Handling the near-planar case
      if (e > EMIN) {
        // Near-planar, elliptical
        double theta = acos(position.x / r);
        double varpi = acos(evec.x / e);
        if (inc < M_PI/2) {
          f = theta - varpi;
        }
        else {
          f = varpi - theta;
        }
      }
      else {
        // Near-planar, near-circular
        f = acos(position.x / r);
        if (velocity.x > 0) {
          f = 2 * M_PI - f;
        }
      }
    } 
    else {
      // Handling the non-planar case
      if (e > EMIN) {
        // Non-planar, elliptical
        double edotr = evec.x*position.x + evec.y*position.y + evec.z*position.z;
        double rdotv = position.x*velocity.x + position.y*velocity.y + position.z*velocity.z;
        f = acos(edotr / (e * r));
        if (rdotv < 0) {
          f = 2 * M_PI - f;
        }
      }
      else {
        // Non-planar, circular
        f = acos((nvec.x*position.x + nvec.y*position.y) / (n * r));
        if (position.z < 0) {
          f = 2 * M_PI - f;
        }
      }
    }
  }
  else {
    double E = acosh((1 - r / a) / e);
    double rdotv = position.x*velocity.x + position.y*velocity.y + position.z*velocity.z;
    if(rdotv < 0) {
     E *= -1;
    }
    // double l = fabs(a * (e*e - 1));
    // f = acos((l / r - 1) / e);

    f = 2 * atan2(sqrt(e + 1) * tanh(E / 2), sqrt(e - 1));
  }

  kep.a = a;
  kep.e = e;
  kep.inc = inc;
  kep.arg = arg;
  kep.node = node;
  kep.f = f;

  return kep;

}

double* py_calc_kep_from_xyz(int N, double mu, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs) {
  //struct KeplerOrbit* output = malloc(N * sizeof(struct KeplerOrbit));
  double* output = malloc(N * 6 * sizeof(double));
  int dummy;
  struct KeplerOrbit kep;

  for (int idx = 0; idx < N; idx++) {
    kep = calc_kep_from_xyz(mu, xs[idx], ys[idx], zs[idx], vxs[idx], vys[idx], vzs[idx]);
    dummy = idx * 6;

    output[dummy] = kep.a;
    output[dummy + 1] = kep.e;
    output[dummy + 2] = kep.inc;
    output[dummy + 3] = kep.arg;
    output[dummy + 4] = kep.node;
    output[dummy + 5] = kep.f;

  }
  return output;
}

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

double* py_calc_M_from_E(int N, double* es, double* Es) {

  double* output = malloc(N * sizeof(double));

  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_M_from_E(es[idx], Es[idx]);
  }

  return output;

}

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

double* py_calc_E_from_f(int N, double* es, double* fs) {

  double* output = malloc(N * sizeof(double));

  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_E_from_f(es[idx], fs[idx]);
  }

  return output;

}


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

double* py_calc_f_from_E(int N, double* es, double* Es) {

  double* output = malloc(N * sizeof(double));

  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_f_from_E(es[idx], Es[idx]);
  }

  return output;

}


double* py_kepM_to_xyz(int N, double *as, double *es, double *incs, double *args, double *nodes, double *Ms)
{

  double* output = malloc(N * 6 * sizeof(double));
  int dummy;

  struct StateVector rock;

  for (int idx = 0; idx < N; idx++) {

    rock = kepM_to_xyz(as[idx], es[idx], incs[idx], args[idx], nodes[idx], Ms[idx]);

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

double* py_kepE_to_xyz(int N, double *as, double *es, double *incs, double *args, double *nodes, double *Es)
{

  double* output = malloc(N * 6 * sizeof(double));
  int dummy;

  struct StateVector rock;

  for (int idx = 0; idx < N; idx++) {

    rock = kepE_to_xyz(as[idx], es[idx], incs[idx], args[idx], nodes[idx], Es[idx]);

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


double* py_calc_E_from_M(int N, double* es, double* Ms) {

  double* output = malloc(N * sizeof(double));

  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_E_from_M(es[idx], Ms[idx]);
  }

  return output;

}

double* py_calc_vovec_from_kep(int N, double mu, double* as, double* es, double* rs, double* Es){

  double* output = malloc(N * 3 * sizeof(double));
  struct Vector3 vec;
  int dummy;

  for (int idx = 0; idx < N; idx++) {
    vec = calc_vovec_from_kep(mu, as[idx], es[idx], rs[idx], Es[idx]);

    dummy = idx * 3;
    output[dummy] = vec.x;
    output[dummy + 1] = vec.y;
    output[dummy + 2] = vec.z;

  }

  return output;

}