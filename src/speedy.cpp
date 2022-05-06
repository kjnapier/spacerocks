#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <thread>
#include <vector>
#include <omp.h>

const double EMIN           = 1e-8;
const double IMIN           = 1e-8;
const double mu_bary        = 0.00029630927493457475;
const double speed_of_light = 173.14463267424034;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// const double ONE_SIXTH = 1.0 / 6.0;

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
        if(fabs(delta_i3) < 1.e-16){
          break;
        }
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

    f = 2 * atan2(sqrt((1 + e)/(1 - e)) * sin(E / 2), cos(E / 2));
    r = a * omece;

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sin(E);
    voy = c * sqrt(1 - e * e) * cosE;

  } else {

    E = calc_E_from_M(e, M);

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

// struct KeplerOrbit calc_kep_from_xyz(double mu, double x, double y, double z, double vx, double vy, double vz) {

//   struct KeplerOrbit kep;
  
//   double r = sqrt(x*x + y*y + z*z);
//   double vsq = vx*vx + vy*vy + vz*vz;
//   double a, e, inc, arg, node, f;
//   double sin_node, cos_node;

//   a = 1 / (2 / r - vsq / mu);

//   double hx = y * vz - z * vy;
//   double hy = z * vx - x * vz;
//   double hz = x * vy - y * vx;

//   double h2 = hx*hx + hy*hy + hz*hz;
//   double h = sqrt(h2);


//   e = sqrt(1 - h2 / (mu * a));
//   inc = acos(hz / h);

//   double sin_i = sin(inc);

//   if (hz > 0) {
//     sin_node = hx / (h * sin_i);
//     cos_node = -hy / (h * sin_i);
//   } else {
//     sin_node = -hx / (h * sin_i);
//     cos_node = hy / (h * sin_i);
//   }
//   node = atan2(sin_node, cos_node);

//   double sin_omega_plus_f = z/(r * sin_i);
//   double cos_omega_plus_f = (x/r + sin_node * cos(inc) * sin_omega_plus_f) / cos_node;

//   double omega_plus_f = atan2(sin_omega_plus_f, cos_omega_plus_f);
//   double p = a * (1 - e*e);

//   double r_dot_rdot = x*vx + y*vy + z*vz;
//   double rdot;
//   if (r_dot_rdot > 0) {
//     rdot = sqrt(vsq + h2/(r*r));
//   } else {
//     rdot = -sqrt(vsq + h2/(r*r));
//   }
  
//   f = atan2(p * rdot / h, p/r - 1);
//   if (f < 0) {
//     f += 2 * M_PI;
//   }

//   arg = omega_plus_f - f;

//   kep.a = a;
//   kep.e = e;
//   kep.inc = inc;
//   kep.arg = arg;
//   kep.node = node;
//   kep.f = f;

//   return kep;

// }


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

  // double rinv = 1 / r;
  // double muinv = 1 / mu;

  // Compute the eccentricity vector
  evec.x = (velocity.y * hvec.z - velocity.z * hvec.y) / mu - position.x / r;
  evec.y = (velocity.z * hvec.x - velocity.x * hvec.z) / mu - position.y / r;
  evec.z = (velocity.x * hvec.y - velocity.y * hvec.x) / mu - position.z / r;

  // evec.x = (velocity.y * hvec.z - velocity.z * hvec.y) * muinv - position.x * rinv;
  // evec.y = (velocity.z * hvec.x - velocity.x * hvec.z) * muinv - position.y * rinv;
  // evec.z = (velocity.x * hvec.y - velocity.y * hvec.x) * muinv - position.z * rinv;


  // Compute the normal vector
  nvec.x = -hvec.y;
  nvec.y = hvec.x;
  double n = sqrt(nvec.x*nvec.x + nvec.y*nvec.y);

  // Compute the semi-major axis (a)
  a = 1 / (2 / r - vsq / mu);
  // a = 1 / (2 * rinv - vsq * muinv);

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

extern "C" {
void free_memory(double *ptr) {
  free(ptr);
}
}

extern "C" {
double* py_calc_kep_from_xyz(int N, double mu, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs) {
  //struct KeplerOrbit* output = malloc(N * sizeof(struct KeplerOrbit));
  double* output = (double*) malloc(N * 6 * sizeof(double));
  //double* output = new double[N * 6];
  int dummy;
  struct KeplerOrbit kep;

  #pragma omp parallel for private(kep, dummy) schedule(guided)
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
}}

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

extern "C" {
double* py_calc_M_from_E(int N, double* es, double* Es) {

  double* output = (double*) malloc(N * sizeof(double));
  //double* output = new double[N];

  #pragma omp parallel for schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_M_from_E(es[idx], Es[idx]);
  }

  return output;

}}

// If inc > pi/2, varpi = node - arg

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

extern "C" {
double* py_calc_E_from_f(int N, double* es, double* fs) {

  double* output = (double*) malloc(N * sizeof(double));
  //double* output = new double[N];

  #pragma omp parallel for default(shared) schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_E_from_f(es[idx], fs[idx]);
  }

  return output;

}}


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

extern "C" {
double* py_calc_f_from_E(int N, double* es, double* Es) {

  double* output = (double*) malloc(N * sizeof(double));
  //double* output = new double[N];

  #pragma omp parallel for default(shared) schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_f_from_E(es[idx], Es[idx]);
  }

  return output;

}}



extern "C" {
double* py_kepM_to_xyz(int N, double *as, double *es, double *incs, double *args, double *nodes, double *Ms)
{

  double* output = (double*) malloc(N * 6 * sizeof(double));
  //double* output = new double[N * 6];
  int dummy;

  struct StateVector rock;

  #pragma omp parallel for private(rock, dummy) schedule(guided)
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

}}

extern "C" {
double* py_kepE_to_xyz(int N, double *as, double *es, double *incs, double *args, double *nodes, double *Es)
{

  double* output = (double*) malloc(N * 6 * sizeof(double));
  //double* output = new double[N * 6];
  int dummy;

  struct StateVector rock;

  #pragma omp parallel for private(rock, dummy) schedule(guided)
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

}}

extern "C" {
double* py_calc_E_from_M(int N, double* es, double* Ms) {

  double* output = (double*) malloc(N * sizeof(double));
  //double* output = new double[N];

  #pragma omp parallel for default(shared) schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    output[idx] = calc_E_from_M(es[idx], Ms[idx]);
  }

  return output;

}}

extern "C" {
double* py_calc_vovec_from_kep(int N, double mu, double* as, double* es, double* rs, double* Es){

  double* output = (double*) malloc(N * 3 * sizeof(double));
  //double* output = new double[N * 3];
  struct Vector3 vec;
  int dummy;

  #pragma omp parallel for private(vec, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    vec = calc_vovec_from_kep(mu, as[idx], es[idx], rs[idx], Es[idx]);

    dummy = idx * 3;
    output[dummy] = vec.x;
    output[dummy + 1] = vec.y;
    output[dummy + 2] = vec.z;

  }

  return output;

}}

// struct StateVector correct_for_ltt(double a, double e, double inc, double arg, double node, double M0, 
//                                   double ox, double oy, double oz, double ovx, double ovy, double ovz) {

//  struct StateVector rock;
//  struct StateVector out;
//  double ltt0 = 0;
//  double dx, dy, dz, dvx, dvy, dvz;
//  double x0, y0, z0;
//  double delta, ltt, dltt;

//  rock = kepM_to_xyz(a, e, inc, arg, node, M0);
//  double n = sqrt(mu_bary / fabs(a*a*a));

//  for (int idx = 0; idx < 5; idx++) {

//    x0 = rock.x;
//    y0 = rock.y;
//    z0 = rock.z;

//    dx = x0 - ox;
//    dy = y0 - oy;
//    dz = z0 - oz;
 
//    delta = sqrt(dx*dx + dy*dy + dz*dz);

//    ltt = delta / speed_of_light;
//    dltt = fabs(ltt - ltt0);

//    if (dltt < 1e-6) {
//      break;
//    }
//    else {
//      rock = kepM_to_xyz(a, e, inc, arg, node, M0 - (ltt * n));
//      ltt0 = ltt;
//    }
//  }

//  dvx = rock.vx - ovx;
//  dvy = rock.vy - ovy;
//  dvz = rock.vz - ovz;

//  out.x  = dx;
//  out.y  = dy;
//  out.z  = dz;
//  out.vx = dvx;
//  out.vy = dvy;
//  out.vz = dvz;

//  return out;

// }

// struct StateVector correct_for_ltt(double a, double e, double inc, double arg, double node, double M0, 
//                                    double ox, double oy, double oz, double ovx, double ovy, double ovz) {

//   struct StateVector rock;
//   struct StateVector temp;
//   struct StateVector out;

//   double ltt0 = 0;
//   double dx, dy, dz, dvx, dvy, dvz;
//   double acc;

//   rock = kepM_to_xyz(a, e, inc, arg, node, M0);
//   double r = sqrt(rock.x*rock.x + rock.y*rock.y + rock.z*rock.z);

//   temp.x  = rock.x;
//   temp.y  = rock.y;
//   temp.z  = rock.z;
//   temp.vx = rock.vx;
//   temp.vy = rock.vy;
//   temp.vz = rock.vz;

//   double xi = mu_bary / (r * r * r);

//   for (int idx = 0; idx < 3; idx++) {

//     dx  = temp.x - ox;
//     dy  = temp.y - oy;
//     dz  = temp.z - oz;
    
//     double delta = sqrt(dx*dx + dy*dy + dz*dz);

//     double ltt = delta / speed_of_light;
//     double dltt = fabs(ltt - ltt0);

//     if (dltt < 1e-6) {
//       break;
//     }
//     else {
      
//       acc = xi * ltt;

//       temp.x = rock.x - (0.5 * acc * rock.x + rock.vx) * ltt;
//       temp.y = rock.y - (0.5 * acc * rock.y + rock.vy) * ltt;
//       temp.z = rock.z - (0.5 * acc * rock.z + rock.vz) * ltt;

//       ltt0 = ltt;
//     }
//   }

//   temp.vx = rock.vx + acc * rock.x;
//   temp.vy = rock.vy + acc * rock.y;
//   temp.vz = rock.vz + acc * rock.z;

//   dvx  = temp.vx - ovx;
//   dvy  = temp.vy - ovy;
//   dvz  = temp.vz - ovz;

//   out.x  = dx;
//   out.y  = dy;
//   out.z  = dz;
//   out.vx = dvx;
//   out.vy = dvy;
//   out.vz = dvz;

//   return out;

// }

struct StateVector correct_for_ltt(double x, double y, double z, double vx, double vy, double vz, 
                                   double ox, double oy, double oz, double ovx, double ovy, double ovz) {

  struct StateVector rock;
  struct StateVector temp;
  struct StateVector out;

  double ltt0 = 0;
  double dx, dy, dz, dvx, dvy, dvz;
  double acc;

  rock.x = x;
  rock.y = y;
  rock.z = z;
  rock.vx = vx;
  rock.vy = vy;
  rock.vz = vz;
  double r = sqrt(rock.x*rock.x + rock.y*rock.y + rock.z*rock.z);

  temp.x  = rock.x;
  temp.y  = rock.y;
  temp.z  = rock.z;
  temp.vx = rock.vx;
  temp.vy = rock.vy;
  temp.vz = rock.vz;

  double xi = mu_bary / (r * r * r);

  for (int idx = 0; idx < 3; idx++) {

    dx  = temp.x - ox;
    dy  = temp.y - oy;
    dz  = temp.z - oz;
    
    double delta = sqrt(dx*dx + dy*dy + dz*dz);

    double ltt = delta / speed_of_light;
    double dltt = fabs(ltt - ltt0);

    if (dltt < 1e-6) {
      break;
    }
    else {
      
      acc = xi * ltt;

      temp.x = rock.x - (0.5 * acc * rock.x + rock.vx) * ltt;
      temp.y = rock.y - (0.5 * acc * rock.y + rock.vy) * ltt;
      temp.z = rock.z - (0.5 * acc * rock.z + rock.vz) * ltt;

      ltt0 = ltt;
    }
  }

  temp.vx = rock.vx + acc * rock.x;
  temp.vy = rock.vy + acc * rock.y;
  temp.vz = rock.vz + acc * rock.z;

  dvx  = temp.vx - ovx;
  dvy  = temp.vy - ovy;
  dvz  = temp.vz - ovz;

  out.x  = dx;
  out.y  = dy;
  out.z  = dz;
  out.vx = dvx;
  out.vy = dvy;
  out.vz = dvz;

  return out;

}

// struct StateVector correct_for_ltt(double x, double y, double z, double vx, double vy, double vz, 
//                                    double ox, double oy, double oz, double ovx, double ovy, double ovz) {

//   struct StateVector temp;
//   struct StateVector out;

//   double ltt0 = 0;
//   double dx, dy, dz, dvx, dvy, dvz;
//   double acc;

//   double r = sqrt(x*x + y*y + z*z);

//   temp.x  = x;
//   temp.y  = y;
//   temp.z  = z;
//   temp.vx = vx;
//   temp.vy = vy;
//   temp.vz = vz;

//   double xi = mu_bary / (r * r * r);

//   for (int idx = 0; idx < 3; idx++) {

//     dx  = temp.x - ox;
//     dy  = temp.y - oy;
//     dz  = temp.z - oz;
    
//     double delta = sqrt(dx*dx + dy*dy + dz*dz);

//     double ltt = delta / speed_of_light;
//     double dltt = fabs(ltt - ltt0);

//     if (dltt < 1e-6) {
//       break;
//     }
//     else {
      
//       acc = xi * ltt;

//       temp.x = x - (0.5 * acc * x + vx) * ltt;
//       temp.y = y - (0.5 * acc * y + vy) * ltt;
//       temp.z = z - (0.5 * acc * z + vz) * ltt;

//       ltt0 = ltt;
//     }
//   }

//   temp.vx = vx + acc * x;
//   temp.vy = vy + acc * y;
//   temp.vz = vz + acc * z;

//   dvx  = temp.vx - ovx;
//   dvy  = temp.vy - ovy;
//   dvz  = temp.vz - ovz;

//   out.x  = dx;
//   out.y  = dy;
//   out.z  = dz;
//   out.vx = dvx;
//   out.vy = dvy;
//   out.vz = dvz;

//   return out;

// }

// extern "C" {
// double* py_correct_for_ltt(int N, double* as, double* es, double* incs, double* args, double* nodes, double* Ms, 
//                            double* obsx, double* obsy, double* obsz, double* obsvx, double* obsvy, double* obsvz) {

//   double* output = (double*) malloc(N * sizeof(struct StateVector));
//   //double* output = new double[N * 6];
//   struct StateVector rock;
//   int dummy;

//   #pragma omp parallel for private(rock, dummy) schedule(guided)
//   for (int idx = 0; idx < N; idx++) {

//     double a    = as[idx];
//     double e    = es[idx];
//     double inc  = incs[idx];
//     double arg  = args[idx];
//     double node = nodes[idx];
//     double M0   = Ms[idx];

//     double ox  = obsx[idx];
//     double oy  = obsy[idx];
//     double oz  = obsz[idx];
//     double ovx = obsvx[idx];
//     double ovy = obsvy[idx];
//     double ovz = obsvz[idx];

//     rock = correct_for_ltt(a, e, inc, arg, node, M0, ox, oy, oz, ovx, ovy, ovz);

//     dummy = idx * 6;
//     output[dummy]     = rock.x;
//     output[dummy + 1] = rock.y;
//     output[dummy + 2] = rock.z;
//     output[dummy + 3] = rock.vx;
//     output[dummy + 4] = rock.vy;
//     output[dummy + 5] = rock.vz;

//   }

//    return output;

//  }}

// extern "C" {
// double* py_correct_for_ltt_single_observer(int N, double* as, double* es, double* incs, double* args, double* nodes, double* Ms, 
//                            double* obsx, double* obsy, double* obsz, double* obsvx, double* obsvy, double* obsvz) {

//   double* output = (double*) malloc(N * sizeof(struct StateVector));
//   //double* output = new double[N * 6];
//   struct StateVector rock;
//   int dummy;

//   double ox  = obsx[0];
//   double oy  = obsy[0];
//   double oz  = obsz[0];
//   double ovx = obsvx[0];
//   double ovy = obsvy[0];
//   double ovz = obsvz[0];

//   #pragma omp parallel for private(rock, dummy) schedule(guided)
//   for (int idx = 0; idx < N; idx++) {

//     double a    = as[idx];
//     double e    = es[idx];
//     double inc  = incs[idx];
//     double arg  = args[idx];
//     double node = nodes[idx];
//     double M0   = Ms[idx];

//     rock = correct_for_ltt(a, e, inc, arg, node, M0, ox, oy, oz, ovx, ovy, ovz);

//     dummy = idx * 6;
//     output[dummy]     = rock.x;
//     output[dummy + 1] = rock.y;
//     output[dummy + 2] = rock.z;
//     output[dummy + 3] = rock.vx;
//     output[dummy + 4] = rock.vy;
//     output[dummy + 5] = rock.vz;

//   }

//    return output;

//  }}

extern "C" {
double* py_correct_for_ltt(int N, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs, 
                           double* obsx, double* obsy, double* obsz, double* obsvx, double* obsvy, double* obsvz) {

  double* output = (double*) malloc(N * sizeof(struct StateVector));
  //double* output = new double[N * 6];
  struct StateVector rock;
  int dummy;

  #pragma omp parallel for private(rock, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {

    double x    = xs[idx];
    double y    = ys[idx];
    double z    = zs[idx];
    double vx   = vxs[idx];
    double vy   = vys[idx];
    double vz   = vzs[idx];

    double ox  = obsx[idx];
    double oy  = obsy[idx];
    double oz  = obsz[idx];
    double ovx = obsvx[idx];
    double ovy = obsvy[idx];
    double ovz = obsvz[idx];

    rock = correct_for_ltt(x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz);

    dummy = idx * 6;
    output[dummy]     = rock.x;
    output[dummy + 1] = rock.y;
    output[dummy + 2] = rock.z;
    output[dummy + 3] = rock.vx;
    output[dummy + 4] = rock.vy;
    output[dummy + 5] = rock.vz;

  }

   return output;

 }}

extern "C" {
double* py_correct_for_ltt_single_observer(int N, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs, 
                                           double* obsx, double* obsy, double* obsz, double* obsvx, double* obsvy, double* obsvz) {

  double* output = (double*) malloc(N * sizeof(struct StateVector));
  //double* output = new double[N * 6];
  struct StateVector rock;
  int dummy;

  double ox  = obsx[0];
  double oy  = obsy[0];
  double oz  = obsz[0];
  double ovx = obsvx[0];
  double ovy = obsvy[0];
  double ovz = obsvz[0];

  #pragma omp parallel for private(rock, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {

    double x    = xs[idx];
    double y    = ys[idx];
    double z    = zs[idx];
    double vx   = vxs[idx];
    double vy   = vys[idx];
    double vz   = vzs[idx];

    rock = correct_for_ltt(x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz);

    dummy = idx * 6;
    output[dummy]     = rock.x;
    output[dummy + 1] = rock.y;
    output[dummy + 2] = rock.z;
    output[dummy + 3] = rock.vx;
    output[dummy + 4] = rock.vy;
    output[dummy + 5] = rock.vz;

  }

   return output;

 }}

extern "C" {
double* py_correct_for_ltt_destnosim(int N, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs, 
                                     double* obsx, double* obsy, double* obsz, double* obsvx, double* obsvy, double* obsvz) {

  double* output = (double*) malloc(N * sizeof(struct Vector3));
  //double* output = new double[N * 6];
  struct StateVector rock;
  int dummy;

  double ox = obsx[0];
  double oy = obsy[0];
  double oz = obsz[0];
  double ovx = obsvx[0];
  double ovy = obsvy[0];
  double ovz = obsvz[0];

  #pragma omp parallel for private(rock, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {

    double x    = xs[idx];
    double y    = ys[idx];
    double z    = zs[idx];
    double vx   = vxs[idx];
    double vy   = vys[idx];
    double vz   = vzs[idx];

    rock = correct_for_ltt(x, y, z, vx, vy, vz, ox, oy, oz, ovx, ovy, ovz);

    dummy = idx * 3;
    output[dummy]     = rock.x;
    output[dummy + 1] = rock.y;
    output[dummy + 2] = rock.z;

  }

   return output;

 }}
