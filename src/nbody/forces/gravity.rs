use crate::nbody::forces::Force;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::state::{position, State};

use nalgebra::Vector3;

/// Implementation of classical Newtonian gravitational force.
///
/// Calculates gravitational interactions between all pairs of bodies using Newton's
/// law of universal gravitation.
#[derive(Debug, Clone, Copy)]
pub struct NewtonianGravity;

impl Force for NewtonianGravity {
    fn add_acceleration(&self, states: &[State], masses: &[f64], acc: &mut [Vector3<f64>]) {
        // Naive implementation of Newtonian gravity. O(0.5 * n^2) complexity. Massive bodies
        // come first, so the loop stops at the first massless one: test particles don't pull
        // on each other.
        for idx in 0..states.len() {
            let m_idx = masses[idx];
            if m_idx == 0.0 {
                break;
            }
            let x_idx = position(&states[idx]);
            let mut acc_idx = Vector3::zeros();
            for jdx in idx + 1..states.len() {
                let r_vec = x_idx - position(&states[jdx]);
                let r2 = r_vec.norm_squared();
                let xi = (GRAVITATIONAL_CONSTANT / (r2 * r2.sqrt())) * r_vec;
                acc_idx -= masses[jdx] * xi;
                acc[jdx] += m_idx * xi;
            }
            acc[idx] += acc_idx;
        }
    }

    fn is_newtonian_gravity(&self) -> bool {
        true
    }
}
