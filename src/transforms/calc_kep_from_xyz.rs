use crate::constants::MU_BARY;
use crate::StateVector;
use crate::KeplerOrbit;
use nalgebra::Vector3;

fn acos_safe(x: f64) -> f64 {
  if x.abs() > 1.0 {
    if x > 0.0 {
      return 0.0;
    }
    else {
      return std::f64::consts::PI;
    }
  }
  else {
    return x.acos();
  }
}

// enum OrbitType {
//   Ellipse,
//   Parabola,
//   Hyperbola,
// }

pub fn calc_kep_from_xyz(state: StateVector) -> KeplerOrbit {

    let EMIN = 1e-10;
    let IMIN = 1e-10;

    let (a, e, inc, mut arg, mut node, mut f);

    let r = state.position.norm();
    let vsq = state.velocity.dot(&state.velocity);

    let hvec = state.position.cross(&state.velocity);
    let evec = state.velocity.cross(&hvec) / MU_BARY - state.position / r;

    let nvec = Vector3::new(-hvec.y, hvec.x, 0.0);
    let n = nvec.norm();

    // let specific_energy = vsq / 2.0 - MU_BARY / r;

    // cases depending on the specific energy
    // < 0: ellipse, 0: parabola, > 0: hyperbola
    // orbit_type = match specific_energy {
    //   e if e < 0.0 => OrbitType::Ellipse,
    //   e if e == 0.0 => OrbitType::Parabola,
    //   e if e > 0.0 => OrbitType::Hyperbola,
    //   _ => panic!("Invalid specific energy"),
    // };

    a = 1.0 / (2.0 / r - vsq / MU_BARY);
    e = evec.norm();
    inc = (hvec.z / hvec.norm()).acos();

    if inc == 0.0 {
      node = 0.0;
    }
    else {
      node = (nvec.x / n).acos();
    }
    if nvec.y < 0.0 {
      node = 2.0 * std::f64::consts::PI - node;
    }
  
    // Compute the argument of pericenter (arg)
    if e < EMIN {
      arg = 0.0;
    }
    else if inc < IMIN || inc > std::f64::consts::PI - IMIN {
      // let mut theta = (state.position.x / r).acos();
      // if state.position.y < 0.0 {
      //   theta = 2.0 * std::f64::consts::PI - theta;
      // }
      let mut varpi = (evec.x / e).acos();
      if evec.y < 0.0 {
        varpi = 2.0 * std::f64::consts::PI - varpi;
      }
      if inc < std::f64::consts::PI / 2.0 {
        arg = varpi - node;
      }
      else {
        arg = node - varpi;
      }
    }
    else {
      let vv = nvec.dot(&evec) / (n * e);
      if (vv.abs() - 1.0).abs() < 1e-10 {
        if vv > 0.0 {
          arg = 0.0;
        }
        else {
          arg = std::f64::consts::PI;
        }
      }
      else {
        arg = vv.acos();
        if evec.z < 0.0 {
          arg = 2.0 * std::f64::consts::PI - arg;
        }
      }      
    }
      
    // Compute the true anomaly (f)
    if e < 1.0 {
      if inc < IMIN || inc > std::f64::consts::PI - IMIN {
        // Handling the near-planar case
        if e > EMIN {
          // Near-planar, elliptical
          let theta = (state.position.x / r).acos();
          let varpi = (evec.x / e).acos();
          if inc < std::f64::consts::PI/2.0 {
            f = theta - varpi;
          }
          else {
            f = varpi - theta;
          }
        }
        else {
          // Near-planar, near-circular
          f = (state.position.x / r).acos();
          if state.velocity.x > 0.0 {
            f = 2.0 * std::f64::consts::PI - f;
          }
        }
      } 
      else {
        // Handling the non-planar case
        if e > EMIN {
          // Non-planar, elliptical
          let edotr = evec.dot(&state.position);
          let rdotv = state.position.dot(&state.velocity);
          let vv = edotr / (e * r);
          if (vv.abs() - 1.0).abs() < 1e-10 {
            if vv > 0.0 {
              f = 0.0;
            }
            else {
              f = std::f64::consts::PI;
            }
          }
          else {
            f = vv.acos();
            if rdotv < 0.0 {
              f = 2.0 * std::f64::consts::PI - f;
            }
          }
          
        }
        else {
          // Non-planar, circular
          f = (nvec.dot(&state.position) / (n * r)).acos();
          if state.position.z < 0.0 {
            f = 2.0 * std::f64::consts::PI - f;
          }
        }
      }
    }
    else {
      let mut argument = (1.0 - r / a) / e;
      if (argument - 1.0).abs() < 1e-10 {
        argument = 1.0;
      }
      let mut E = argument.acosh();
      let rdotv = state.position.dot(&state.velocity);
      if rdotv < 0.0 {
       E *= -1.0;
      }
      f = 2.0 * ((e + 1.0).sqrt() * (E / 2.0).tanh()).atan2((e - 1.0).sqrt());
    }
  
    let kep = KeplerOrbit::new(a, e, inc, arg, node, f);
    return kep;
      
}