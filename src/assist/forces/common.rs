//! Helpers shared by the force models.

use nalgebra::{Matrix3, Vector3};

use crate::assist::spice_simulation::{SimulationParticle, SimulationState};
use crate::constants::GRAVITATIONAL_CONSTANT;

/// NAIF codes of the Sun, planets (barycenters), Earth, Moon and Pluto, in ASSIST's order
/// (the bodies ASSIST calls `ASSIST_BODY_NPLANETS`).
pub const MAJOR_BODIES: [i32; 11] = [10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9];

pub const SUN: i32 = 10;
pub const EARTH: i32 = 399;

/// Index of the perturber with NAIF code `code` in `state.spice_bodies`.
#[inline]
pub fn body_index(state: &SimulationState, code: i32) -> Option<usize> {
    state.spice_bodies.iter().position(|b| b.code == code)
}

/// GM (AU^3/day^2), position and velocity of perturber `i`.
#[inline]
pub fn body(state: &SimulationState, i: usize) -> (f64, Vector3<f64>, Vector3<f64>) {
    let p = &state.spice_particles[i];
    (GRAVITATIONAL_CONSTANT * p.mass, p.position, p.velocity)
}

/// Add position and velocity partials to a particle's acceleration Jacobian.
#[inline]
pub fn add_jacobian(p: &mut SimulationParticle, dadr: &Matrix3<f64>, dadv: Option<&Matrix3<f64>>) {
    for i in 0..3 {
        for j in 0..3 {
            p.stm[3 + i][j] += dadr[(i, j)];
        }
    }
    if let Some(dv) = dadv {
        for i in 0..3 {
            for j in 0..3 {
                p.stm[3 + i][3 + j] += dv[(i, j)];
            }
        }
    }
}

/// Rotation from J2000 to a body's equatorial frame with pole (right ascension, declination),
/// in the convention ASSIST uses (x' = pole x z_J2000 direction, z' = pole).
pub fn pole_rotation(ra: f64, dec: f64) -> Matrix3<f64> {
    let (sina, cosa) = ra.sin_cos();
    let (sind, cosd) = dec.sin_cos();
    Matrix3::new(
        -sina, cosa, 0.0,
        -cosa * sind, -sina * sind, cosd,
        cosa * cosd, sina * cosd, sind,
    )
}
