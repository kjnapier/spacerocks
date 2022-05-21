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

double calc_E_from_M(double e, double M);
double calc_E_from_f(double e, double f);
double calc_M_from_E(double e, double E);
double calc_f_from_E(double e, double E);

struct Vector3 calc_vovec_from_kep(double mu, double a, double e, double r, double E);
struct StateVector kepM_to_xyz(double a, double e, double inc, double arg, double node, double M);
struct StateVector kepE_to_xyz(double a, double e, double inc, double arg, double node, double E);
struct StateVector kepf_to_xyz(double a, double e, double inc, double arg, double node, double f);
struct KeplerOrbit calc_kep_from_xyz(double mu, double x, double y, double z, double vx, double vy, double vz);

struct StateVector correct_for_ltt(double x, double y, double z, double vx, double vy, double vz, 
                                   double ox, double oy, double oz, double ovx, double ovy, double ovz);

