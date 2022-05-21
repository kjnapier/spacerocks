#include "spacerocks.hpp"

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