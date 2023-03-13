#include "spacerocks.hpp"

extern "C" {
void free_memory(double *ptr) {
  free(ptr);
}}

extern "C" {
double* py_compute_topocentric_correction(int N, double* lats, double* lons, double* elevations, double* epochs) {

  double* output = (double*) malloc(N * 3 * sizeof(double));
  struct Vector3 vec;
  int dummy;

  #pragma omp parallel for private(vec, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    vec = compute_topocentric_correction(lats[idx], lons[idx], elevations[idx], epochs[idx]);
    dummy = idx * 3;

    output[dummy] = vec.x;
    output[dummy + 1] = vec.y;
    output[dummy + 2] = vec.z;

  }
  return output;
}}


// extern "C" {
// double* py_calc_kep_from_xyz(int N, double mu, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs) {
//   //struct KeplerOrbit* output = malloc(N * sizeof(struct KeplerOrbit));
//   double* output = (double*) malloc(N * 6 * sizeof(double));
//   //double* output = new double[N * 6];
//   int dummy;
//   struct KeplerOrbit kep;

//   #pragma omp parallel for private(kep, dummy) schedule(guided)
//   for (int idx = 0; idx < N; idx++) {
//     kep = calc_kep_from_xyz(mu, xs[idx], ys[idx], zs[idx], vxs[idx], vys[idx], vzs[idx]);
//     dummy = idx * 6;

//     output[dummy] = kep.a;
//     output[dummy + 1] = kep.e;
//     output[dummy + 2] = kep.inc;
//     output[dummy + 3] = kep.arg;
//     output[dummy + 4] = kep.node;
//     output[dummy + 5] = kep.f;

//   }
//   return output;
// }}

extern "C" {
void py_calc_kep_from_xyz(int N, double mu, double* xs, double* ys, double* zs, double* vxs, double* vys, double* vzs, double* output) {

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
  return;
}}


// extern "C" {
// double* py_calc_M_from_E(int N, double* es, double* Es) {

//   double* output = (double*) malloc(N * sizeof(double));
//   //double* output = new double[N];

//   #pragma omp parallel for schedule(guided)
//   for (int idx = 0; idx < N; idx++) {
//     output[idx] = calc_M_from_E(es[idx], Es[idx]);
//   }

//   return output;

// }}

extern "C" {
void py_calc_M_from_E(int N, double* es, double* Es, double* Ms) {

  #pragma omp parallel for schedule(guided)
  for (int idx = 0; idx < N; idx++) {
    Ms[idx] = calc_M_from_E(es[idx], Es[idx]);
  }
  return;
}}

// If inc > pi/2, varpi = node - arg

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
double* py_kepf_to_xyz(int N, double *as, double *es, double *incs, double *args, double *nodes, double *fs)
{

  double* output = (double*) malloc(N * 6 * sizeof(double));
  //double* output = new double[N * 6];
  int dummy;

  struct StateVector rock;

  #pragma omp parallel for private(rock, dummy) schedule(guided)
  for (int idx = 0; idx < N; idx++) {

    rock = kepf_to_xyz(as[idx], es[idx], incs[idx], args[idx], nodes[idx], fs[idx]);

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
