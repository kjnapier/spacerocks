use crate::nbody::forces::{central_body, Force};
use crate::constants::{GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT};
use crate::state::{pv, State};

use nalgebra::Vector3;

#[derive(Debug, Clone, Copy)]
/// Post-Newtonian correction to the gravity of the central body (the most massive one, the
/// Sun in a solar system simulation).
pub struct SolarGR;

impl Force for SolarGR {
    fn add_acceleration(&self, states: &[State], masses: &[f64], acc: &mut [Vector3<f64>]) {
        let Some(sun) = central_body(masses) else { return };
        let mu = GRAVITATIONAL_CONSTANT * masses[sun];
        let (sun_x, sun_v) = pv(&states[sun]);

        for (idx, s) in states.iter().enumerate() {
            if idx == sun {
                continue;
            }
            let (x, v) = pv(s);
            let r_vec = x - sun_x;
            let r = r_vec.norm();

            let v_vec = v - sun_v;
            let v = v_vec.norm();

            let s0 = mu / (SPEED_OF_LIGHT.powi(2) * r * r * r);
            let s1 = ((4.0 * mu) / r - v * v) * r_vec;
            let s2 = 4.0 * (r_vec.dot(&v_vec)) * v_vec;

            acc[idx] += s0 * (s1 + s2);
        }
    }
}
