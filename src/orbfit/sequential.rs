//! Keeping orbits up to date as detections arrive: layup's sequential (information-filter)
//! update and its incremental driver (`sequential_update`, `incremental_orbitfit`; issue #419).
//!
//! A converged fit summarizes the detections it used by its parameters and covariance. When new
//! detections arrive, [`sequential_update`] fits only them, with that summary as a Gaussian prior
//! (its information matrix added to the normal equations). As long as the prior is close to
//! Gaussian this equals a refit of everything, at the cost of integrating only the new
//! detections. A large move (in prior sigmas) means the linearization can't be trusted, and the
//! update falls back to a full refit.
//!
//! [`update_orbit`] routes one object the way layup's driver does: skip it if its detections are
//! unchanged, update it sequentially if detections were only added, and refit it otherwise.

use nalgebra::DMatrix;

use crate::orbfit::astrometry::fingerprint;
use crate::orbfit::lm::{fit_orbit, fit_orbit_with_prior};
use crate::orbfit::pipeline::determine_orbit;
use crate::orbfit::{Astrometry, FitFlag, FitOptions, OrbitFit};
use crate::spice::SpiceKernel;

/// How [`update_orbit`] handled an object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UpdateRoute {
    /// The detections are the ones the prior was fitted to: the prior is returned as it is.
    Skip,
    /// Detections were only added, and the sequential update was accepted.
    Sequential,
    /// Detections were only added, but the sequential update failed or moved too far, so every
    /// detection was refitted from the prior.
    SequentialFallback,
    /// Detections were removed or changed: every detection was refitted from the prior.
    Full,
    /// No usable prior: fitted from scratch.
    Cold,
}

impl UpdateRoute {
    /// layup's routing names: "skip", "sequential", "sequential_fallback", "full", "cold".
    pub fn name(&self) -> &'static str {
        match self {
            UpdateRoute::Skip => "skip",
            UpdateRoute::Sequential => "sequential",
            UpdateRoute::SequentialFallback => "sequential_fallback",
            UpdateRoute::Full => "full",
            UpdateRoute::Cold => "cold",
        }
    }
}

/// The fitted parameters of `fit` (state, then its fitted non-gravitational parameters).
fn parameters(fit: &OrbitFit) -> Vec<f64> {
    let mut p = fit.state.to_vec();
    p.extend((0..3).filter(|&k| fit.fit_nongrav[k]).map(|k| fit.nongrav[k]));
    p
}

fn covariance(fit: &OrbitFit) -> DMatrix<f64> {
    DMatrix::from_row_slice(fit.npar, fit.npar, &fit.covariance)
}

/// Mahalanobis distance of `updated` from `prior` in the prior's sigmas (layup's
/// `_update_mahalanobis`, over every fitted parameter). Infinite if the covariance is singular.
pub fn update_mahalanobis(prior: &OrbitFit, updated: &OrbitFit) -> f64 {
    let (p0, p1) = (parameters(prior), parameters(updated));
    if p0.len() != p1.len() || p0.len() != prior.npar {
        return f64::INFINITY;
    }
    let dx = nalgebra::DVector::from_iterator(p0.len(), p1.iter().zip(&p0).map(|(a, b)| a - b));
    match covariance(prior).lu().solve(&dx) {
        Some(s) => dx.dot(&s).max(0.0).sqrt(),
        None => f64::INFINITY,
    }
}

/// Update the converged fit `prior` with the detections `new` alone (layup's
/// `sequential_update`). The epoch stays the prior's, and the parameters the prior fitted (its
/// state, and any non-gravitational parameters) are fitted again, starting from the prior, with
/// the prior's covariance as a Gaussian prior ([`fit_orbit_with_prior`]).
///
/// If that fit fails (including a prior covariance that isn't positive definite, flag 7) or moves
/// the parameters more than `opts.max_update_sigma` prior sigmas, and `all` (every detection, old
/// and new) is given, the result is instead a fit of `all` starting from the prior. Without `all`
/// the failed update is returned, with flag 8 for too large a move.
///
/// Returns the fit and whether the sequential update was accepted. An accepted update reports
/// the chi-square of `new` alone, with `ndof` its number of residual rows (layup's convention),
/// and the posterior covariance.
pub fn sequential_update(new: &Astrometry, prior: &OrbitFit, all: Option<&Astrometry>, kernel: &SpiceKernel, opts: &FitOptions) -> (OrbitFit, bool) {
    let opts = FitOptions { fit_nongrav: prior.fit_nongrav, nongrav_auto: None, nongrav_per_arc: false, ..opts.clone() };
    let refit = || all.map(|a| fit_orbit(a, prior.epoch, &prior.state, &prior.nongrav, kernel, &opts, opts.max_iter));

    let info = if prior.per_arc || prior.npar != parameters(prior).len() {
        None
    } else {
        covariance(prior).try_inverse().map(|m| (&m + m.transpose()) * 0.5)
    };
    let seq = match info {
        // layup's run_sequential_update: the absolute convergence test
        Some(info) => fit_orbit_with_prior(new, prior.epoch, &prior.state, &prior.nongrav, &info, kernel, &FitOptions { conv_frac: 0.0, ..opts.clone() }, opts.max_iter),
        None => OrbitFit { epoch: prior.epoch, state: prior.state, ..OrbitFit::failed(FitFlag::PriorNotPositiveDefinite) },
    };
    if seq.flag != FitFlag::Converged {
        return match refit() {
            Some(f) => (f, false),
            None => (seq, false),
        };
    }
    if update_mahalanobis(prior, &seq) > opts.max_update_sigma {
        return match refit() {
            Some(f) => (f, false),
            None => (OrbitFit { flag: FitFlag::NonlinearUpdate, ..seq }, false),
        };
    }
    (seq, true)
}

/// A previous fit and the keys ([`Astrometry::detection_key`]) of the detections it was fitted
/// to, in any order.
#[derive(Debug, Clone)]
pub struct PriorFit {
    pub fit: OrbitFit,
    pub keys: Vec<u64>,
}

/// Bring one object's orbit up to date with its current detections, as layup's
/// `incremental_orbitfit` routes each object:
///
/// - **skip**: the detections are exactly those of a converged prior (same fingerprint and
///   count). The prior is returned, its residuals and `used` put in the current order.
/// - **sequential**: detections were only added. [`sequential_update`] with the new ones, falling
///   back to a refit of all of them.
/// - **full**: a converged prior exists, but detections were removed or changed.
///   [`determine_orbit`] of all detections, starting from the prior.
/// - **cold**: no converged prior. [`determine_orbit`] from scratch.
///
/// After a sequential update the residuals are filled in for the new detections only (NaN for
/// the others), since the old ones were not integrated.
pub fn update_orbit(current: &Astrometry, prior: Option<&PriorFit>, kernel: &SpiceKernel, opts: &FitOptions) -> (OrbitFit, UpdateRoute) {
    let keys = current.detection_keys();
    let prior = prior.filter(|p| p.fit.converged());
    let Some(p) = prior else {
        return (determine_orbit(current, None, kernel, opts), UpdateRoute::Cold);
    };

    let n = current.len();
    if p.keys.len() == n && fingerprint(&p.keys) == fingerprint(&keys) {
        let mut fit = p.fit.clone();
        let mut by_key: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
        for (j, &k) in p.keys.iter().enumerate() {
            by_key.entry(k).or_default().push(j);
        }
        let order: Vec<usize> = keys.iter().map(|k| by_key.get_mut(k).and_then(|v| v.pop()).unwrap_or(0)).collect();
        if p.fit.residuals.len() == n {
            fit.residuals = order.iter().map(|&j| p.fit.residuals[j]).collect();
        }
        fit.used = if p.fit.used.len() == n { order.iter().map(|&j| p.fit.used[j]).collect() } else { vec![true; n] };
        return (fit, UpdateRoute::Skip);
    }

    let old: std::collections::HashSet<u64> = p.keys.iter().copied().collect();
    let now: std::collections::HashSet<u64> = keys.iter().copied().collect();
    if old.is_subset(&now) {
        let new_idx: Vec<usize> = (0..n).filter(|&i| !old.contains(&keys[i])).collect();
        if !new_idx.is_empty() {
            let (mut fit, accepted) = sequential_update(&current.subset(&new_idx), &p.fit, Some(current), kernel, opts);
            if accepted {
                let mut residuals = vec![[f64::NAN; 6]; n];
                for (r, &i) in fit.residuals.iter().zip(&new_idx) {
                    residuals[i] = *r;
                }
                fit.residuals = residuals;
                fit.used = vec![true; n];
                return (fit, UpdateRoute::Sequential);
            }
            return (fit, UpdateRoute::SequentialFallback);
        }
    }
    let fit = determine_orbit(current, Some((p.fit.epoch, p.fit.state, p.fit.nongrav)), kernel, opts);
    (fit, UpdateRoute::Full)
}
