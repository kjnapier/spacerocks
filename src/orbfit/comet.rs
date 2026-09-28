//! Original and future orbits of long-period comets (layup's `comet`).
//!
//! The "original" orbit is the barycentric osculating orbit a comet had on its way in, before
//! the planets perturbed it, and the "future" one the orbit it leaves on. Both are conventionally
//! evaluated at 250 AU from the barycenter (as in the CODE catalogue and layup), where the planets
//! act as a point mass at the barycenter and the barycentric two-body elements (with the mass of
//! the Sun and planets) no longer change.
//!
//! The orbit is integrated with the full force model to the inbound (original) or outbound
//! (future) crossing of the reference distance. If the planetary ephemeris runs out first (a
//! comet reaches 250 AU centuries from perihelion), the elements are taken where the
//! integration stopped, which is exact in the two-body approximation that holds out there; layup
//! instead continues with the planets as REBOUND particles. `distance` says where the elements
//! were evaluated and `reached` whether that was the reference distance.

use nalgebra::Vector3;

use crate::orbfit::residuals::trial_simulation;
use crate::orbfit::{FitOptions, MU_TOTAL};
use crate::spice::SpiceKernel;

/// layup's reference distance for original and future orbits (AU).
pub const REFERENCE_DISTANCE: f64 = 250.0;

/// Barycentric osculating two-body elements (GM of the Sun and planets) at some point.
#[derive(Debug, Clone, Copy)]
pub struct CometOrbit {
    /// TDB Julian date where the elements were evaluated.
    pub epoch: f64,
    /// Barycentric distance there (AU).
    pub distance: f64,
    /// Whether that is the reference distance (false: the ephemeris ran out first).
    pub reached: bool,
    /// 1/a (1/AU; negative for hyperbolic orbits), a (AU), e, q (AU), inclination (radians,
    /// to the J2000 ecliptic).
    pub inv_a: f64,
    pub a: f64,
    pub e: f64,
    pub q: f64,
    pub inc: f64,
}

/// Barycentric two-body elements of `state` (J2000 equatorial) with GM `mu`.
pub fn barycentric_elements(state: &[f64; 6], mu: f64) -> (f64, f64, f64, f64) {
    let r = Vector3::new(state[0], state[1], state[2]);
    let v = Vector3::new(state[3], state[4], state[5]);
    let inv_a = 2.0 / r.norm() - v.norm_squared() / mu;
    let h = r.cross(&v);
    let e_vec = v.cross(&h) / mu - r / r.norm();
    let e = e_vec.norm();
    let q = h.norm_squared() / mu / (1.0 + e);
    // inclination to the J2000 ecliptic
    let eps = 84381.448_f64.to_radians() / 3600.0;
    let hz = h.z * eps.cos() - h.y * eps.sin();
    let inc = (hz / h.norm()).clamp(-1.0, 1.0).acos();
    (inv_a, e, q, inc)
}

/// Where the comet is (barycentric distance and radial velocity) and its elements.
fn orbit_here(t: f64, state: &[f64; 6], reached: bool) -> CometOrbit {
    let (inv_a, e, q, inc) = barycentric_elements(state, MU_TOTAL);
    let r = (state[0] * state[0] + state[1] * state[1] + state[2] * state[2]).sqrt();
    CometOrbit { epoch: t, distance: r, reached, inv_a, a: 1.0 / inv_a, e, q, inc }
}

/// The original (`future = false`) or future orbit of the comet with barycentric J2000 `state`
/// at TDB `epoch` and non-gravitational parameters `nongrav` (with `opts.gofr`), at barycentric
/// distance `reference` (AU). `None` if the orbit never gets that far: a barycentric aphelion
/// inside it, or a perihelion outside it.
pub fn comet_orbit(epoch: f64, state: &[f64; 6], nongrav: &[f64; 3], future: bool, reference: f64, kernel: &SpiceKernel, opts: &FitOptions) -> Option<CometOrbit> {
    let (inv_a, e, q, _) = barycentric_elements(state, MU_TOTAL);
    if q > reference || (inv_a > 0.0 && e < 1.0 && q * (1.0 + e) / (1.0 - e) < reference) {
        return None;
    }
    let radial = |s: &[f64; 6]| (s[0] * s[3] + s[1] * s[4] + s[2] * s[5]) / (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    let dist = |s: &[f64; 6]| (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    // The crossing wanted has rdot < 0 (original) or > 0 (future). Outside the reference distance
    // and already on that leg, it is ahead in time for the original orbit (behind for the future
    // one); otherwise it is the other way.
    let (r0, v0) = (dist(state), radial(state));
    let on_leg = if future { v0 > 0.0 } else { v0 < 0.0 };
    let forward = if r0 >= reference && on_leg { !future } else { future };
    if r0 >= reference && on_leg && (r0 - reference).abs() < 1e-9 {
        return Some(orbit_here(epoch, state, true));
    }

    let mut sim = trial_simulation(epoch, state, nongrav, &[], kernel, opts, false).ok()?;
    let jd_ref = sim.state.jd_ref;
    let sign = if forward { 1.0 } else { -1.0 };
    let get = |sim: &crate::assist::SpiceSimulation| {
        let p = &sim.state.particles[0];
        [p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z]
    };
    let on_wanted_leg = |s: &[f64; 6]| if future { radial(s) > 0.0 } else { radial(s) < 0.0 };
    let mut t = epoch - jd_ref;
    let mut last = *state;
    // March in chunks sized to cover about a tenth of the remaining distance, until r - reference
    // changes sign on the wanted leg; then bisect in time for the crossing.
    loop {
        let speed = (last[3] * last[3] + last[4] * last[4] + last[5] * last[5]).sqrt().max(1e-6);
        let chunk = ((reference - dist(&last)).abs().max(1.0) * 0.1 / speed).clamp(1.0, 5000.0);
        let t_next = t + sign * chunk;
        if sim.integrate_rel(t_next, kernel).is_err() {
            // Out of the ephemeris: the barycentric two-body elements where it stopped.
            return Some(orbit_here(jd_ref + t, &last, false));
        }
        let now = get(&sim);
        let (f0, f1) = (dist(&last) - reference, dist(&now) - reference);
        if f0 * f1 <= 0.0 && on_wanted_leg(&now) {
            let (mut lo, mut hi) = (t, t_next);
            let mut best = now;
            for _ in 0..50 {
                let mid = 0.5 * (lo + hi);
                if sim.integrate_rel(mid, kernel).is_err() {
                    break;
                }
                let sm = get(&sim);
                best = sm;
                if (dist(&sm) - reference) * f0 > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
                if (hi - lo).abs() < 1e-6 {
                    break;
                }
            }
            return Some(orbit_here(jd_ref + 0.5 * (lo + hi), &best, true));
        }
        t = t_next;
        last = now;
        if t.abs() > 400_000.0 {
            return Some(orbit_here(jd_ref + t, &last, false));
        }
    }
}

/// Both [`comet_orbit`]s at [`REFERENCE_DISTANCE`]: (original, future).
pub fn original_and_future(epoch: f64, state: &[f64; 6], nongrav: &[f64; 3], kernel: &SpiceKernel, opts: &FitOptions) -> (Option<CometOrbit>, Option<CometOrbit>) {
    (
        comet_orbit(epoch, state, nongrav, false, REFERENCE_DISTANCE, kernel, opts),
        comet_orbit(epoch, state, nongrav, true, REFERENCE_DISTANCE, kernel, opts),
    )
}
