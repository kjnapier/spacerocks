use crate::constants::MU_BARY;
use crate::StateVector;
use crate::time::Time;

use crate::transforms::calc_kep_from_xyz;
use crate::transforms::calc_E_from_f;
use crate::transforms::calc_f_from_E;
use crate::transforms::calc_M_from_E;

pub enum UniversalOrbitType {
    Hyperbolic,
    Parabolic,
    Elliptical,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniversalOrbit {
    pub e: f64,
    pub q: f64,
    pub inc: f64,
    pub node: f64,
    pub arg: f64,
    pub universal_anomaly: f64,
    pub mu: f64,
}

impl UniversalOrbit {

    pub fn new(e: f64, q: f64, inc: f64, node: f64, arg: f64, universal_anomaly: f64, mu: f64) -> Result<UniversalOrbit, String> {

        if e < 0.0 {
            return Err("Eccentricity must be greater than or equal to zero.".to_string());
        }

        if q < 0.0 {
            return Err("Pericenter distance must be greater than or equal to zero.".to_string());
        }

        UniversalOrbit {
            e: e,
            q: q,
            inc: inc,
            node: node,
            arg: arg,
            universal_anomaly: universal_anomaly,
            mu: mu,
        }

    }


}
