use crate::constants::MU_BARY;
use crate::StateVector;

use crate::transforms::calc_kep_from_xyz;
use crate::transforms::calc_E_from_f;
use crate::transforms::calc_f_from_E;
use crate::transforms::calc_M_from_E;

use serde::{Serialize, Deserialize};

// TODO: This should require mu as an argument

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KeplerOrbit {
    pub a: f64,
    pub e: f64,
    pub inc: f64,
    pub arg: f64,
    pub node: f64,
    pub f: f64,
}

impl KeplerOrbit {
    pub fn from_xyz(state: StateVector) -> Self {
        let kep = calc_kep_from_xyz::calc_kep_from_xyz(state);
        KeplerOrbit {
            a: kep.a,
            e: kep.e,
            inc: kep.inc,
            arg: kep.arg,
            node: kep.node,
            f: kep.f,
        }
    }

    pub fn new(a: f64, e: f64, inc: f64, arg: f64, node: f64, f: f64) -> Self {
        KeplerOrbit {
            a: a,
            e: e,
            inc: inc,
            arg: arg,
            node: node,
            f: f,
        }
    }

    pub fn n(&self) -> f64 {
        (MU_BARY / (self.a.powi(3)).abs()).sqrt()
    }

    pub fn M(&self) -> f64 {
        let E = self.E();
        let M = calc_M_from_E::calc_M_from_E(self.e, E);
        return M;
    }

    pub fn E(&self) -> f64 {
        let E = calc_E_from_f::calc_E_from_f(self.e, self.f);
        return E;
    }

    pub fn varpi(&self) -> f64 {
        return (self.arg + self.node) % (2.0 * std::f64::consts::PI);
    }

    pub fn q(&self) -> f64 {
        return self.a * (1.0 - self.e);
    }

}
