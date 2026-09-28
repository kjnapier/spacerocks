//! Orbit determination and fitting from optical astrometry.
//!
//! A port of the orbit fitter in layup (Smithsonian/layup): the same initial orbit
//! determination, candidate screening, Levenberg–Marquardt iteration and force model (ASSIST's,
//! via [`crate::assist::SpiceSimulation`]), so that a fit here reproduces layup's to within the
//! integrator's tolerance.
//!
//! The data are an [`Astrometry`] (parallel arrays of epochs, RA, Dec, uncertainties and
//! barycentric observer positions); orbits are barycentric J2000 states `[f64; 6]` (AU, AU/day)
//! at a TDB Julian date.
//!
//! - [`determine_orbit`]: the whole pipeline (IOD, candidate selection, differential
//!   correction), or a differential correction from a given initial orbit.
//! - [`fit_orbit`]: the differential correction alone.
//! - [`residuals`]: residuals and partial derivatives of a trial orbit.
//! - [`gauss_states`] / [`gauss`]: Gauss's method.
//! - [`bk_iod`]: Bernstein–Khushalani linear initial orbit, the fallback when Gauss fails.
//! - [`herget_iod`]: Herget's method, iterating the ranges to the first and last detections.
//! - [`veres_sigma`]: astrometric uncertainties after Vereš et al. (2017), as layup assigns them.
//! - [`BiasTable`]: star-catalog debiasing (Eggl et al. 2020), as layup applies it.

pub mod astrometry;
    pub use astrometry::{ades_observer_state, fingerprint, geodetic_to_earth_fixed, occultation_radec, observer_accelerations, observer_positions, observer_states, Astrometry, RowKind, ARCSEC, DEFAULT_DELAY_SIGMA, DEFAULT_RATE_SIGMA, DEFAULT_SIGMA};

pub mod gauss;
    pub use gauss::{gauss, gauss_astrometry, gauss_states, MU_TOTAL};

pub mod residuals;
    pub use residuals::{residuals, residuals_per_arc, Residuals};

pub mod lm;
    pub use lm::{fit_orbit, fit_orbit_per_arc, fit_orbit_with_prior};

pub mod bk;
    pub use bk::bk_iod;

pub mod bk_fit;
    pub use bk_fit::fit_orbit_bk;

pub mod herget;
    pub use herget::{herget, herget_iod};

pub mod weights;
    pub use weights::veres_sigma;

pub mod debias;
    pub use debias::BiasTable;

pub mod predict;
    pub use predict::{predict, predict_astrometry, Prediction};

pub mod comet;
    pub use comet::{comet_orbit, original_and_future, CometOrbit};

pub mod sequential;
    pub use sequential::{sequential_update, update_mahalanobis, update_orbit, PriorFit, UpdateRoute};

pub mod pipeline;
    pub use pipeline::{build_sequence, determine_orbit, select_nongrav};

/// Outcome of a fit. The values are layup's `flag` codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitFlag {
    /// Not attempted: fewer than three detections, a non-finite value, or a detection before
    /// 1801.
    NotAttempted = -1,
    /// Converged, and every check passed.
    Converged = 0,
    /// The differential correction did not converge.
    NotConverged = 1,
    /// Reduced chi-square above [`FitOptions::chi2_threshold`].
    Chi2TooLarge = 2,
    /// Initial orbits were found, but none converged on the primary arc.
    NoRootConverged = 3,
    /// The primary arc converged, but the incremental extension to all detections did not.
    BuildupFailed = 4,
    /// No initial orbit could be determined (neither Gauss's method nor the
    /// Bernstein–Khushalani fallback produced a candidate).
    NoSolution = 5,
    /// Converged, but a fitted non-gravitational parameter has no usable variance.
    DegenerateCovariance = 6,
    /// Sequential update: the prior covariance is not positive definite (and no full set of
    /// detections was given to refit instead).
    PriorNotPositiveDefinite = 7,
    /// Sequential update: the new detections moved the orbit more than
    /// `FitOptions::max_update_sigma` prior sigmas (and no full set was given to refit instead).
    NonlinearUpdate = 8,
    /// Converged, but the hyperbolic excess speed exceeds 200 km/s.
    ImplausibleOrbit = 9,
}

impl FitFlag {
    pub fn code(&self) -> i32 {
        *self as i32
    }

    /// The flag with layup's code `code`, if there is one.
    pub fn from_code(code: i32) -> Option<FitFlag> {
        use FitFlag::*;
        [NotAttempted, Converged, NotConverged, Chi2TooLarge, NoRootConverged, BuildupFailed, NoSolution, DegenerateCovariance, PriorNotPositiveDefinite, NonlinearUpdate, ImplausibleOrbit]
            .into_iter()
            .find(|f| f.code() == code)
    }

    pub fn description(&self) -> &'static str {
        match self {
            FitFlag::NotAttempted => "not attempted (fewer than 3 valid detections)",
            FitFlag::Converged => "converged",
            FitFlag::NotConverged => "did not converge",
            FitFlag::Chi2TooLarge => "reduced chi-square too large",
            FitFlag::NoRootConverged => "no initial orbit converged",
            FitFlag::BuildupFailed => "incremental fit to all detections failed",
            FitFlag::NoSolution => "no initial orbit",
            FitFlag::DegenerateCovariance => "non-gravitational parameter unconstrained",
            FitFlag::PriorNotPositiveDefinite => "prior covariance not positive definite",
            FitFlag::NonlinearUpdate => "update too large for the prior's linearization",
            FitFlag::ImplausibleOrbit => "implausible hyperbolic excess speed",
        }
    }
}

/// Thresholds of the automatic non-gravitational model selection (layup's
/// `NongravAutoThresholds`, for `fit_nongrav="auto"`). The defaults are layup's.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NongravAuto {
    /// A gravity-only fit with chi-square per degree of freedom at or below this is kept, and
    /// no non-gravitational model is tried. 0 always tries them.
    pub accept_reduced_chi2: f64,
    /// A model must lower chi-square by more than this per added parameter.
    pub delta_chi2_per_param: f64,
    /// Every added parameter must exceed this many times its 1-sigma uncertainty.
    pub nsigma: f64,
}

impl Default for NongravAuto {
    fn default() -> Self {
        NongravAuto { accept_reduced_chi2: 1.5, delta_chi2_per_param: 9.0, nsigma: 3.0 }
    }
}

/// Settings for the fit. The defaults are layup's.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// Iteration limit of a differential correction.
    pub max_iter: usize,
    /// Iteration limit when screening initial-orbit candidates.
    pub screen_iter: usize,
    /// Convergence: every component of the last step below this (AU, AU/day).
    pub tolerance: f64,
    /// Scaled convergence (layup's `conv_frac`, issue #477): with a positive value, a parameter
    /// has converged when its step is below `max(tolerance, conv_frac * sigma)`, sigma being its
    /// formal uncertainty. 0 (the default) keeps the absolute test.
    pub conv_frac: f64,
    /// Largest acceptable chi-square per degree of freedom.
    pub chi2_threshold: f64,
    /// IAS15 tolerance.
    pub epsilon: f64,
    /// Which non-gravitational parameters (A1, A2, A3) to fit. They are fitted after the
    /// gravity-only orbit has converged, starting from it.
    pub fit_nongrav: [bool; 3],
    /// Piecewise-constant non-gravitational parameters (layup's `per_arc`, for linking comet
    /// apparitions): the detections before the fit epoch (arc A) and after it (arc B) each get
    /// their own values of the parameters in `fit_nongrav`, with one shared state. The fit epoch
    /// must lie between the two apparitions. Not used with `nongrav_auto`.
    pub nongrav_per_arc: bool,
    /// Choose the non-gravitational model automatically (layup's `fit_nongrav="auto"`); overrides
    /// `fit_nongrav`. See [`pipeline::select_nongrav`].
    pub nongrav_auto: Option<NongravAuto>,
    /// Marsden g(r) as `[alpha, m, n, k, r0]`, `g = alpha (r/r0)^-m (1 + (r/r0)^n)^-k`.
    /// `None` is the inverse-square law used for asteroids.
    pub gofr: Option<[f64; 5]>,
    /// Detections separated by more than this many days start a new arc.
    pub arc_gap: f64,
    /// Initial orbits closer than this to the barycenter (AU) are used only if nothing else
    /// converges.
    pub min_distance: f64,
    /// Initial orbits whose 80th-percentile residual exceeds this many sigma are discarded
    /// before fitting.
    pub prefilter_sigma: f64,
    /// When no Gauss candidate converges, try a Bernstein–Khushalani linear seed (layup's
    /// `iod="auto"`; `false` is its `iod="gauss"`).
    pub bk_fallback: bool,
    /// Seed from Herget's method on the primary arc instead of Gauss's (layup's
    /// `iod="herget"`). There is no Bernstein–Khushalani fallback on this path.
    pub herget: bool,
    /// Run the forward and backward integrations of a residual evaluation in parallel.
    pub parallel: bool,
    /// Robust mode (not in layup; off by default). Rejects outliers, and if layup's pipeline
    /// finds no orbit, starts over from the short window of detections that gives one and
    /// widens it. See [`determine_orbit`].
    pub robust: bool,
    /// Robust mode: detections whose normalized residual (rms over their measurements)
    /// exceeds this many sigma are left out of the fit.
    pub outlier_sigma: f64,
    /// Robust mode: the longest window (days) an initial orbit is sought in.
    pub seed_window: f64,
    /// Robust mode: how many windows (most detections first) to try.
    pub seed_tries: usize,
    /// The fitting engine for the pipeline's gravity-only fits (candidate screening, the fit to
    /// all detections, and the arc-by-arc build-up): layup's `engine`.
    pub engine: Engine,
    /// Sequential updates: the largest move, in prior sigmas (Mahalanobis), accepted before
    /// falling back to a full refit (layup's `max_update_sigma`).
    pub max_update_sigma: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            max_iter: 100,
            screen_iter: 80,
            tolerance: 1e-12,
            conv_frac: 0.0,
            chi2_threshold: 10.0,
            epsilon: 1e-9,
            fit_nongrav: [false; 3],
            nongrav_auto: None,
            nongrav_per_arc: false,
            gofr: None,
            arc_gap: 90.0,
            min_distance: 0.3,
            prefilter_sigma: 1000.0,
            bk_fallback: true,
            herget: false,
            parallel: true,
            robust: false,
            outlier_sigma: 4.0,
            seed_window: 60.0,
            seed_tries: 10,
            max_update_sigma: 4.0,
            engine: Engine::Cartesian,
        }
    }
}

/// layup's fitting engines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Engine {
    /// Levenberg–Marquardt in barycentric Cartesian coordinates ([`fit_orbit`]).
    #[default]
    Cartesian,
    /// Levenberg–Marquardt in Bernstein–Khushalani parameters with a bound-orbit prior on the
    /// line-of-sight velocity ([`fit_orbit_bk`]). As in layup, fits with non-gravitational
    /// parameters stay Cartesian. It also applies to the fit from an `initial` orbit (layup's
    /// `initial_guess` path is always Cartesian).
    BkNative,
}

/// A fitted orbit.
#[derive(Debug, Clone)]
pub struct OrbitFit {
    /// TDB Julian date of `state`.
    pub epoch: f64,
    /// Barycentric J2000 position (AU) and velocity (AU/day).
    pub state: [f64; 6],
    /// Non-gravitational parameters A1, A2, A3 (AU/day^2).
    pub nongrav: [f64; 3],
    /// Which of `nongrav` were fitted (the others were held fixed).
    pub fit_nongrav: [bool; 3],
    /// With per-arc non-gravitational parameters (`per_arc`), those of arc B, the detections
    /// after the epoch; `nongrav` then holds arc A's.
    pub nongrav_arc2: [f64; 3],
    pub per_arc: bool,
    /// Row-major `npar x npar` covariance of the state and the fitted non-gravitational
    /// parameters (in that order; with `per_arc`, arc A's then arc B's): the inverse of the
    /// normal matrix at the last iteration.
    pub covariance: Vec<f64>,
    pub npar: usize,
    pub chi2: f64,
    pub ndof: i64,
    pub niter: usize,
    pub flag: FitFlag,
    /// Residuals at the last iteration, per detection: `[ra, dec, ra_rate, dec_rate, delay,
    /// doppler]`, observed minus computed, NaN where not measured. RA (times cos Dec) and Dec in
    /// radians, rates in radians/day, delay in days, Doppler in AU/day.
    pub residuals: Vec<[f64; 6]>,
    /// Which detections the fit used: all of them, unless robust mode rejected some. Empty for
    /// a fit that produced no orbit.
    pub used: Vec<bool>,
}

impl OrbitFit {
    /// A failed fit with no orbit.
    pub fn failed(flag: FitFlag) -> OrbitFit {
        OrbitFit {
            epoch: f64::NAN,
            state: [f64::NAN; 6],
            nongrav: [0.0; 3],
            fit_nongrav: [false; 3],
            nongrav_arc2: [0.0; 3],
            per_arc: false,
            covariance: vec![f64::NAN; 36],
            npar: 6,
            chi2: f64::NAN,
            ndof: 0,
            niter: 0,
            flag,
            residuals: Vec::new(),
            used: Vec::new(),
        }
    }

    pub fn converged(&self) -> bool {
        self.flag == FitFlag::Converged
    }

    /// The 6x6 covariance of the state.
    pub fn state_covariance(&self) -> [[f64; 6]; 6] {
        let mut c = [[f64::NAN; 6]; 6];
        if self.covariance.len() == self.npar * self.npar {
            for (i, row) in c.iter_mut().enumerate() {
                for (j, v) in row.iter_mut().enumerate() {
                    *v = self.covariance[i * self.npar + j];
                }
            }
        }
        c
    }

    /// 1-sigma uncertainties of the fitted non-gravitational parameters (NaN for those held
    /// fixed; arc A's with `per_arc`).
    pub fn nongrav_sigma(&self) -> [f64; 3] {
        self.amplitude_sigma(6)
    }

    /// With `per_arc`, the 1-sigma uncertainties of arc B's non-gravitational parameters.
    pub fn nongrav_arc2_sigma(&self) -> [f64; 3] {
        if !self.per_arc {
            return [f64::NAN; 3];
        }
        self.amplitude_sigma(6 + self.fit_nongrav.iter().filter(|&&b| b).count())
    }

    fn amplitude_sigma(&self, first: usize) -> [f64; 3] {
        let mut s = [f64::NAN; 3];
        let mut j = first;
        for k in 0..3 {
            if self.fit_nongrav[k] && j < self.npar && self.covariance.len() == self.npar * self.npar {
                s[k] = self.covariance[j * self.npar + j].sqrt();
                j += 1;
            }
        }
        s
    }
}
