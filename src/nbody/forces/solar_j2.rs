use crate::nbody::forces::{central_body, Force};
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::state::{position, State};

use nalgebra::Vector3;

/// Implementation of perturbations due to the Sun's oblateness (J2).
/// Accounts for the Sun's slight equatorial bulge which creates
/// a non-spherical component to its gravitational field. The Sun is the central body (the
/// most massive one).
#[derive(Debug, Clone, Copy)]
pub struct SolarJ2;

const SUN_J2: f64 = 2.17e-7;
const SUN_RADIUS: f64 = 696_342.0 / 149_597_870.7;

impl Force for SolarJ2 {
    /// Acceleration due to the Sun's oblateness (J₂). This perturbation depends on latitude
    /// with respect to the solar equator (through the z-component) and falls off as r⁻⁵.
    fn add_acceleration(&self, states: &[State], masses: &[f64], acc: &mut [Vector3<f64>]) {
        let Some(sun) = central_body(masses) else { return };
        let mu = GRAVITATIONAL_CONSTANT * masses[sun];
        let sun_x = position(&states[sun]);

        for (idx, s) in states.iter().enumerate() {
            if idx == sun {
                continue;
            }
            let r_vec = position(s) - sun_x;
            let r = r_vec.norm();
            let z2_r2 = r_vec[2] * r_vec[2] / (r * r);

            let factor = 2.0 * SUN_J2 * mu * SUN_RADIUS * SUN_RADIUS / (2.0 * r.powi(5));
            let ax = factor * r_vec.x * (5.0 * z2_r2 - 1.0);
            let ay = factor * r_vec.y * (5.0 * z2_r2 - 1.0);
            let az = factor * r_vec.z * (5.0 * z2_r2 - 3.0);

            acc[idx] += Vector3::new(ax, ay, az);
        }
    }
}
