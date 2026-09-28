//! Predicted positions of a fitted orbit, with uncertainty (layup's `predict`).
//!
//! At each time the orbit is integrated to the light-time-corrected emission time for the
//! observer, together with its variational equations. The partials of the line of sight with
//! respect to the fitted parameters (including the light-time term, as layup's) are projected on
//! the local sky basis (A along increasing RA, D along increasing Dec) and map the fit's
//! covariance to a 2x2 on-sky covariance: `B C B^T`. Unlike layup, which maps the 6x6 state
//! covariance and propagates gravity-only uncertainty, the covariance here is the fit's full one,
//! so fitted non-gravitational parameters (and per-arc ones) add their uncertainty too.

use nalgebra::Vector3;

use crate::orbfit::astrometry::Astrometry;
use crate::orbfit::residuals::{integrate_light_time, trial_simulation};
use crate::orbfit::{FitOptions, OrbitFit};
use crate::spice::SpiceKernel;

const SPEED_OF_LIGHT: f64 = 173.1446326742403;

/// One predicted position.
#[derive(Debug, Clone, Copy)]
pub struct Prediction {
    /// TDB Julian date (receive time).
    pub epoch: f64,
    /// Astrometric (light-time corrected) RA and Dec, radians.
    pub ra: f64,
    pub dec: f64,
    /// Distance from the observer at the emission time, AU.
    pub delta: f64,
    /// On-sky covariance (radians²) in the (RA·cos Dec, Dec) directions.
    pub cov: [[f64; 2]; 2],
}

impl Prediction {
    /// The 1-sigma error ellipse: semi-major and semi-minor axes (radians) and the position
    /// angle of the major axis (degrees, North through East, in [0, 180)), as layup's
    /// `skyplane_cov_to_radec_cov`.
    pub fn ellipse(&self) -> (f64, f64, f64) {
        let [[xx, xy], [_, yy]] = self.cov;
        let tr = xx + yy;
        let disc = ((xx - yy).powi(2) + (2.0 * xy).powi(2)).sqrt();
        let a = (0.5 * (tr + disc)).max(0.0).sqrt();
        let b = (0.5 * (tr - disc)).max(0.0).sqrt();
        let pa = (90.0 - 0.5 * (2.0 * xy).atan2(xx - yy).to_degrees()).rem_euclid(180.0);
        (a, b, pa)
    }
}

/// Predict `fit` at the TDB Julian dates `epochs` for observers at the barycentric J2000
/// positions `observers` (AU). The orbit keeps its fitted non-gravitational parameters (per arc
/// if it was fitted that way).
pub fn predict(fit: &OrbitFit, epochs: &[f64], observers: &[[f64; 3]], kernel: &SpiceKernel, opts: &FitOptions) -> Result<Vec<Prediction>, String> {
    Ok(predict_with_states(fit, epochs, observers, kernel, opts, true)?.into_iter().map(|(p, _)| p).collect())
}

/// [`predict`], also returning the object's barycentric J2000 state at each emission time.
/// Without `partials` no variational equations are integrated and the covariances are zero.
pub(crate) fn predict_with_states(
    fit: &OrbitFit,
    epochs: &[f64],
    observers: &[[f64; 3]],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<Vec<(Prediction, [f64; 6])>, String> {
    if epochs.len() != observers.len() {
        return Err(format!("{} epochs for {} observers", epochs.len(), observers.len()));
    }
    if !fit.state.iter().all(|v| v.is_finite()) {
        return Err("the fit has no orbit".into());
    }
    let active: Vec<usize> = (0..3).filter(|&k| fit.fit_nongrav[k]).collect();
    let nact = active.len();
    let npar = fit.npar;
    let cov = |i: usize, j: usize| fit.covariance.get(i * npar + j).copied().unwrap_or(f64::NAN);

    // Forward (after the epoch, in time order) and backward passes, as the fit's.
    let mut order: Vec<usize> = (0..epochs.len()).collect();
    order.sort_by(|&a, &b| epochs[a].total_cmp(&epochs[b]));
    let forward: Vec<usize> = order.iter().copied().filter(|&i| epochs[i] > fit.epoch).collect();
    let backward: Vec<usize> = order.iter().rev().copied().filter(|&i| epochs[i] <= fit.epoch).collect();

    let opts = FitOptions { fit_nongrav: fit.fit_nongrav, ..opts.clone() };
    let mut out = vec![None; epochs.len()];
    for (seq, arc_b) in [(&forward, true), (&backward, false)] {
        if seq.is_empty() {
            continue;
        }
        let ng = if fit.per_arc && arc_b { fit.nongrav_arc2 } else { fit.nongrav };
        // Where this pass's amplitude columns sit in the covariance.
        let amp_at = if fit.per_arc && arc_b { 6 + nact } else { 6 };
        let mut sim = trial_simulation(fit.epoch, &fit.state, &ng, &active, kernel, &opts, partials).map_err(|e| e.to_string())?;
        let jd_ref = sim.state.jd_ref;
        for &i in seq.iter() {
            let r_obs = Vector3::from(observers[i]);
            integrate_light_time(&mut sim, epochs[i] - jd_ref, &r_obs, kernel).map_err(|e| e.to_string())?;
            let p = &sim.state.particles[0];
            let d = p.position - r_obs;
            let dist = d.norm();
            let rho = d / dist;
            let (ra, dec) = (rho.y.atan2(rho.x).rem_euclid(std::f64::consts::TAU), rho.z.asin());
            let a_vec = Vector3::new(-ra.sin(), ra.cos(), 0.0);
            let d_vec = Vector3::new(-dec.sin() * ra.cos(), -dec.sin() * ra.sin(), dec.cos());
            let obj = [p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z];
            if !partials {
                out[i] = Some((Prediction { epoch: epochs[i], ra, dec, delta: dist, cov: [[0.0; 2]; 2] }, obj));
                continue;
            }
            // B: 2 x npar, the sky partials with respect to the covariance's columns.
            let mut b = vec![[0.0f64; 2]; npar];
            for (j, vp) in sim.state.variational_particles.iter().take(6 + nact).enumerate() {
                let ddist = rho.dot(&vp.position);
                let drho = vp.position / dist - rho * (ddist / dist) - p.velocity * (ddist / (dist * SPEED_OF_LIGHT));
                let col = if j < 6 { j } else { amp_at + j - 6 };
                if col < npar {
                    b[col] = [drho.dot(&a_vec), drho.dot(&d_vec)];
                }
            }
            let mut c = [[0.0; 2]; 2];
            for (u, cu) in c.iter_mut().enumerate() {
                for (v, cuv) in cu.iter_mut().enumerate() {
                    *cuv = (0..npar).map(|p1| (0..npar).map(|p2| b[p1][u] * cov(p1, p2) * b[p2][v]).sum::<f64>()).sum();
                }
            }
            out[i] = Some((Prediction { epoch: epochs[i], ra, dec, delta: dist, cov: c }, obj));
        }
    }
    Ok(out.into_iter().map(|p| p.expect("every epoch is in one pass")).collect())
}

/// [`predict`] at the detections of `astrometry` (their epochs and observers).
pub fn predict_astrometry(fit: &OrbitFit, astrometry: &Astrometry, kernel: &SpiceKernel, opts: &FitOptions) -> Result<Vec<Prediction>, String> {
    predict(fit, &astrometry.epoch, &astrometry.observer, kernel, opts)
}
