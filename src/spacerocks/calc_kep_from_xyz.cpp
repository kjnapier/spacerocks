#include "spacerocks.hpp"

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
    double argument = (1 - r / a) / e;
    if (abs(argument - 1) < 1e-10) {
      argument = 1;
    }
    double E = acosh(argument);
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