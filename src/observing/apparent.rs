//! Lightweight observable quantities, computed without allocating.

use nalgebra::Vector3;

use crate::spacerock::hg_magnitude;
use crate::transforms::correct_for_ltt_vectors;

/// Observable quantities of an object seen by an observer, corrected for light travel time.
///
/// This is the plain-number counterpart of [`crate::Observation`] (no `Time`, `Observer` or
/// covariance attached), used by the array and batch APIs. Angles are in radians, rates in
/// radians/day, distances in AU and range rate in AU/day. `mag` is NaN when the object has no
/// absolute magnitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Apparent {
    pub ra: f64,
    pub dec: f64,
    pub ra_rate: f64,
    pub dec_rate: f64,
    /// Observer-object distance at emission time (AU).
    pub range: f64,
    pub range_rate: f64,
    /// Sun-object distance at emission time (AU).
    pub r_helio: f64,
    /// Sun-object-observer angle (radians).
    pub phase: f64,
    /// Sun-observer-object angle (radians).
    pub elong: f64,
    pub mag: f64,
}

impl Apparent {
    /// All fields NaN (used as a fill value).
    pub const NAN: Apparent = Apparent {
        ra: f64::NAN,
        dec: f64::NAN,
        ra_rate: f64::NAN,
        dec_rate: f64::NAN,
        range: f64::NAN,
        range_rate: f64::NAN,
        r_helio: f64::NAN,
        phase: f64::NAN,
        elong: f64::NAN,
        mag: f64::NAN,
    };
}

/// Compute [`Apparent`] quantities from state vectors.
///
/// `position`/`velocity` are the object's state and `obs_position`/`obs_velocity`/`sun_position`
/// the observer's and the Sun's positions, all in the same frame and origin (the light-time
/// correction assumes the origin is the solar system barycenter or the Sun). RA/Dec are
/// measured in that frame's equator, so pass J2000 states for ICRF RA/Dec.
#[inline]
pub fn apparent(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    obs_position: &Vector3<f64>,
    obs_velocity: &Vector3<f64>,
    sun_position: &Vector3<f64>,
    h: Option<f64>,
    g: f64,
) -> Apparent {
    let (p, v) = correct_for_ltt_vectors(position, velocity, obs_position, obs_velocity);

    let mut ra = p.y.atan2(p.x);
    if ra < 0.0 {
        ra += 2.0 * std::f64::consts::PI;
    }
    let dec = (p.z / p.norm()).asin();

    let xi = p.x.powi(2) + p.y.powi(2);
    let ra_rate = -(p.y * v.x - p.x * v.y) / xi;
    let num = -p.z * (p.x * v.x + p.y * v.y) + xi * v.z;
    let denom = xi.sqrt() * p.norm_squared();
    let dec_rate = num / denom;

    let range = p.norm();
    let range_rate = p.dot(&v) / range;

    // Object at emission time relative to the Sun.
    let obj_helio = p + obs_position - sun_position;
    let r_helio = obj_helio.norm();
    let phase = obj_helio.angle(&p);
    let elong = (sun_position - obs_position).angle(&p);

    let mag = match h {
        Some(h) => hg_magnitude(h, g, r_helio, range, phase),
        None => f64::NAN,
    };

    Apparent { ra, dec, ra_rate, dec_rate, range, range_rate, r_helio, phase, elong, mag }
}
