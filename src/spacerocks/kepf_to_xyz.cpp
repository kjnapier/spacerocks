#include "spacerocks.hpp"

struct StateVector kepf_to_xyz(double a, double e, double inc, double arg, double node, double f) {

  double r, c, ox, oy, vox, voy;
  double si, sa, sn, ci, ca, cn;
  double c1, c2, c3, c4, c5, c6;

  struct StateVector rock;

  double E = calc_E_from_f(e, f);

  if (e < 1) {

    
    double cosE = cos(E);
    double omece = 1 - e * cosE;
    r = a * omece;

    c = sqrt(mu_bary * a) / r;

    ox = r * cos(f);
    oy = r * sin(f);
    vox = - c * sin(E);
    voy = c * sqrt(1 - e * e) * cosE;

  } else {

    
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