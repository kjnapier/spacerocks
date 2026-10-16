//! Markov chain Monte Carlo sampling of an orbit's posterior.
//!
//! [`mcmc`] runs a random-walk Metropolis chain over the parameters of a fit (the barycentric
//! state at its epoch, then the fitted non-gravitational parameters), with a flat prior and the
//! likelihood `exp(-chi2 / 2)` of the astrometry under the fit's force model. Proposals are
//! Gaussian, shaped by the fit's covariance and scaled by `2.38 / sqrt(npar)` (Gelman, Roberts
//! & Gilks 1996); during the burn-in the scale adapts toward a 25% acceptance rate. The chain
//! starts `spread` sigmas from the fit (along a random direction), so the burn-in shows.
//!
//! The burn-in is only a view: [`McmcChain::posterior`] drops the first `burn_in` steps, which
//! can be changed after the chain has run.
use nalgebra::{DMatrix, DVector};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::orbfit::{residuals, Astrometry, FitOptions, OrbitFit};
use crate::spice::SpiceKernel;

/// Settings of a chain.
#[derive(Debug, Clone, PartialEq)]
pub struct McmcOptions {
    /// Steps to take (each one residual evaluation).
    pub steps: usize,
    /// Steps to discard from the start when summarizing the posterior.
    pub burn_in: usize,
    /// How far from the fit the chain starts, in sigmas of the fit's covariance.
    pub spread: f64,
    /// Random seed (`None` for a random one).
    pub seed: Option<u64>,
}

impl Default for McmcOptions {
    fn default() -> Self {
        McmcOptions { steps: 5000, burn_in: 1000, spread: 3.0, seed: None }
    }
}

/// One step of a chain.
#[derive(Debug, Clone, PartialEq)]
pub struct McmcStep {
    /// The chain's position after the step: the state, then the fitted non-gravitational
    /// parameters.
    pub params: Vec<f64>,
    pub chi2: f64,
    /// Whether the step's proposal was accepted.
    pub accepted: bool,
}

/// A chain over the parameters of `fit`.
#[derive(Debug, Clone, Default)]
pub struct McmcChain {
    /// TDB Julian date of the states.
    pub epoch: f64,
    /// Which non-gravitational parameters (A1, A2, A3) follow the state in each step's params.
    pub fit_nongrav: [bool; 3],
    /// The values of the non-gravitational parameters held fixed.
    pub nongrav: [f64; 3],
    pub steps: Vec<McmcStep>,
    /// Whether the chain was stopped before taking every step.
    pub stopped: bool,
}

impl McmcChain {
    /// Number of parameters per step.
    pub fn npar(&self) -> usize {
        6 + self.fit_nongrav.iter().filter(|&&b| b).count()
    }

    /// The state of step `i`.
    pub fn state(&self, i: usize) -> [f64; 6] {
        std::array::from_fn(|k| self.steps[i].params[k])
    }

    /// The non-gravitational parameters of step `i` (fitted ones from the chain, others fixed).
    pub fn nongrav_at(&self, i: usize) -> [f64; 3] {
        let mut ng = self.nongrav;
        let mut j = 6;
        for k in 0..3 {
            if self.fit_nongrav[k] {
                ng[k] = self.steps[i].params[j];
                j += 1;
            }
        }
        ng
    }

    /// The steps after the first `burn_in`.
    pub fn posterior(&self, burn_in: usize) -> &[McmcStep] {
        &self.steps[burn_in.min(self.steps.len())..]
    }

    /// Fraction of proposals accepted after the first `burn_in` steps.
    pub fn acceptance(&self, burn_in: usize) -> f64 {
        let post = self.posterior(burn_in);
        post.iter().filter(|s| s.accepted).count() as f64 / post.len().max(1) as f64
    }

    /// Mean and covariance (row-major) of the parameters after the first `burn_in` steps.
    pub fn moments(&self, burn_in: usize) -> (Vec<f64>, Vec<f64>) {
        let post = self.posterior(burn_in);
        let (n, d) = (post.len() as f64, self.npar());
        let mut mean = vec![0.0; d];
        for s in post {
            for k in 0..d {
                mean[k] += s.params[k] / n;
            }
        }
        let mut cov = vec![0.0; d * d];
        for s in post {
            for i in 0..d {
                for j in 0..d {
                    cov[i * d + j] += (s.params[i] - mean[i]) * (s.params[j] - mean[j]) / (n - 1.0).max(1.0);
                }
            }
        }
        (mean, cov)
    }
}

/// A standard normal deviate (Box–Muller).
fn normal(rng: &mut StdRng) -> f64 {
    let u: f64 = 1.0 - rng.gen::<f64>();
    let v: f64 = rng.gen();
    (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}

/// Chi-square of `astrometry` at `state` and `nongrav` (NaN if the orbit can't be integrated).
fn chi2(astrometry: &Astrometry, epoch: f64, state: &[f64; 6], nongrav: &[f64; 3], kernel: &SpiceKernel, opts: &FitOptions) -> f64 {
    match residuals(astrometry, epoch, state, nongrav, kernel, opts, false) {
        Ok(r) => r.resid.iter().zip(r.row_detection.iter().zip(&r.row_kind)).map(|(v, (&i, &k))| (v / astrometry.sigma(i, k)).powi(2)).sum(),
        Err(_) => f64::NAN,
    }
}

/// Sample the posterior of `fit`'s orbit given `astrometry` (the detections it used) with a
/// random-walk Metropolis chain (see the module docs). `opts` gives the force model; the
/// fitted non-gravitational parameters are `fit.fit_nongrav`. `progress` is called after every
/// step with the step number and the chain so far; when it returns false the chain stops.
pub fn mcmc(
    astrometry: &Astrometry,
    fit: &OrbitFit,
    kernel: &SpiceKernel,
    opts: &FitOptions,
    mc: &McmcOptions,
    mut progress: impl FnMut(usize, &McmcChain) -> bool,
) -> Result<McmcChain, String> {
    if !(fit.state.iter().all(|v| v.is_finite()) && fit.epoch.is_finite()) {
        return Err("the fit has no orbit to sample around".into());
    }
    if fit.per_arc {
        return Err("per-arc non-gravitational parameters can't be sampled".into());
    }
    let d = 6 + fit.fit_nongrav.iter().filter(|&&b| b).count();
    if fit.npar != d || fit.covariance.len() != d * d || fit.covariance.iter().any(|v| !v.is_finite()) {
        return Err("the fit has no usable covariance".into());
    }
    let opts = FitOptions { fit_nongrav: fit.fit_nongrav, nongrav_auto: None, nongrav_per_arc: false, progress: None, ..opts.clone() };

    // Proposal shape: the Cholesky factor of the covariance, with a little jitter if needed.
    let cov = DMatrix::from_row_slice(d, d, &fit.covariance);
    let mut jitter = 0.0;
    let l = loop {
        let mut c = cov.clone();
        for i in 0..d {
            c[(i, i)] += jitter * cov[(i, i)].abs().max(f64::MIN_POSITIVE);
        }
        if let Some(ch) = c.cholesky() {
            break ch.l();
        }
        jitter = if jitter == 0.0 { 1e-12 } else { jitter * 100.0 };
        if jitter > 1e-2 {
            return Err("the fit's covariance isn't positive definite".into());
        }
    };

    let mut rng = match mc.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    let gauss = |rng: &mut StdRng| DVector::from_fn(d, |_, _| normal(rng));

    let mut chain = McmcChain { epoch: fit.epoch, fit_nongrav: fit.fit_nongrav, nongrav: fit.nongrav, steps: Vec::with_capacity(mc.steps), stopped: false };
    let mut best = Vec::with_capacity(d);
    best.extend_from_slice(&fit.state);
    for k in 0..3 {
        if fit.fit_nongrav[k] {
            best.push(fit.nongrav[k]);
        }
    }
    let unpack = |p: &DVector<f64>| -> ([f64; 6], [f64; 3]) {
        let state = std::array::from_fn(|k| p[k]);
        let mut ng = fit.nongrav;
        let mut j = 6;
        for k in 0..3 {
            if fit.fit_nongrav[k] {
                ng[k] = p[j];
                j += 1;
            }
        }
        (state, ng)
    };
    let eval = |p: &DVector<f64>| {
        let (s, ng) = unpack(p);
        chi2(astrometry, fit.epoch, &s, &ng, kernel, &opts)
    };

    // Start `spread` sigmas from the fit along a random direction (or at the fit).
    let best = DVector::from_vec(best);
    let mut x = best.clone();
    let mut chi = f64::NAN;
    for attempt in 0..20 {
        let spread = mc.spread / (1 << attempt.min(10)) as f64;
        let z = gauss(&mut rng);
        let norm = z.norm().max(f64::MIN_POSITIVE);
        x = if spread > 0.0 { &best + &l * (z * (spread / norm)) } else { best.clone() };
        chi = eval(&x);
        if chi.is_finite() {
            break;
        }
    }
    if !chi.is_finite() {
        return Err("the orbit can't be integrated near the fit".into());
    }

    let mut scale = 2.38 / (d as f64).sqrt();
    let mut window = 0usize;
    for step in 0..mc.steps {
        let proposal = &x + &l * (gauss(&mut rng) * scale);
        let c = eval(&proposal);
        let accepted = c.is_finite() && (c <= chi || rng.gen::<f64>() < (-(c - chi) / 2.0).exp());
        if accepted {
            x = proposal;
            chi = c;
            window += 1;
        }
        // Adapt the scale toward 25% acceptance, every 50 steps of the burn-in.
        if step < mc.burn_in && (step + 1) % 50 == 0 {
            let rate = window as f64 / 50.0;
            scale *= ((rate - 0.25) * 2.0).exp().clamp(0.5, 2.0);
            window = 0;
        }
        chain.steps.push(McmcStep { params: x.iter().copied().collect(), chi2: chi, accepted });
        if !progress(step, &chain) {
            chain.stopped = true;
            break;
        }
    }
    Ok(chain)
}
