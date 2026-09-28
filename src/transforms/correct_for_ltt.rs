use crate::constants::*;
use crate::StateVector;
use crate::SpaceRock;
use crate::Observer;
use nalgebra::Vector3;

/// Calculates the observer-centric state vector of a rock, accounting for light-time travel.
///
/// # Arguments
/// * `rock` - A SpaceRock object representing the rock.
/// * `observer` - An Observer object representing the observer.
///
/// # Returns
/// * A StateVector object representing the observer-centric state vector of the rock.
///
/// The rock's position at the emission time is found by a second-order Taylor expansion
/// (velocity plus a central-force acceleration), iterated twice on the light time; its velocity
/// at emission is corrected to first order. If the observer has no velocity, it is taken as zero
/// (the position is unaffected; the relative velocity then refers to a stationary observer).
pub fn correct_for_ltt(rock: &SpaceRock, observer: &Observer) -> StateVector {
    let obs_vel = observer.velocity.unwrap_or_else(Vector3::zeros);
    let (d_pos, d_vel) = correct_for_ltt_vectors(&rock.position, &rock.velocity, &observer.position, &obs_vel);
    StateVector::new(d_pos, d_vel)
}

/// Vector form of [`correct_for_ltt`]: returns the observer-centric position and velocity of an
/// object with barycentric state (`r0`, `v0`) seen from an observer at (`obs_pos`, `obs_vel`).
#[inline]
pub fn correct_for_ltt_vectors(
    r0: &Vector3<f64>,
    v0: &Vector3<f64>,
    obs_pos: &Vector3<f64>,
    obs_vel: &Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let r0 = *r0;
    let v0 = *v0;
    let obs_pos = *obs_pos;

    let xi = MU_BARY / r0.norm().powi(3);
    let inv_c = 1.0 / SPEED_OF_LIGHT;

    let mut d_pos = r0 - obs_pos;
    let mut ltt = d_pos.norm() * inv_c;

    // iteration 1
    let mut acc = xi * ltt;
    let pos1 = r0 - v0 * ltt - r0 * (0.5 * acc * ltt);
    d_pos = pos1 - obs_pos;
    ltt = d_pos.norm() * inv_c;

    // iteration 2
    acc = xi * ltt;
    let pos2 = r0 - v0 * ltt - r0 * (0.5 * acc * ltt);
    d_pos = pos2 - obs_pos;
    ltt = d_pos.norm() * inv_c;

    acc = xi * ltt;
    let d_vel = (v0 + r0 * acc) - obs_vel;

    (d_pos, d_vel)
}
