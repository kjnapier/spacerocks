//! Levenberg–Marquardt differential correction, following layup's `orbit_fit`.

use nalgebra::{DMatrix, DVector};

use crate::orbfit::residuals::{residuals, residuals_per_arc};
use crate::orbfit::{Astrometry, FitFlag, FitOptions, OrbitFit};
use crate::spice::SpiceKernel;

/// Seed for fitted non-gravitational parameters that start at zero: the force (and so the
/// variational equations for the parameter) is skipped when all three are zero.
const NONGRAV_SEED: f64 = 1e-20;

/// Arcseconds per radian, as layup rounds it (sets the initial damping).
const ARCSEC_PER_RAD: f64 = 206265.0;

/// Refine the orbit `state` (barycentric J2000, AU and AU/day) at TDB Julian date `epoch` by
/// least squares against `astrometry`, with at most `max_iter` iterations.
///
/// `nongrav` holds the starting (A1, A2, A3); those selected in `opts.fit_nongrav` are fitted
/// along with the state and the others are held fixed.
///
/// The iteration is layup's: each step solves the damped problem
/// `[sqrt(W) B; sqrt(lambda) I] dx = [-sqrt(W) r; 0]` by QR, is accepted when the gain ratio
/// exceeds 0.1 (halving lambda) and otherwise doubles lambda, and the fit has converged when
/// every component of the step is below `opts.tolerance` (or, with `opts.conv_frac > 0` and once
/// a step has been accepted, below `max(tolerance, conv_frac * sigma_i)`, where `sigma_i` is the
/// parameter's formal uncertainty from the current normal matrix: layup's issue #477). A
/// converged fit whose reduced
/// chi-square exceeds `opts.chi2_threshold` is flagged [`FitFlag::Chi2TooLarge`].
///
/// With `opts.nongrav_per_arc` (and parameters to fit), arc A (before `epoch`) and arc B (after)
/// each get their own parameters, both starting from `nongrav`: see [`fit_orbit_per_arc`].
pub fn fit_orbit(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    max_iter: usize,
) -> OrbitFit {
    let arc2 = if opts.nongrav_per_arc && opts.fit_nongrav.iter().any(|&b| b) { Some(*nongrav) } else { None };
    fit(astrometry, epoch, state, nongrav, arc2, None, kernel, opts, max_iter)
}

/// [`fit_orbit`] with a Gaussian prior centred on the starting parameters, given by its
/// information matrix `prior_info` (the inverse covariance; `npar x npar`, state first, then
/// the parameters selected in `opts.fit_nongrav`). This is layup's sequential-update fit (issue
/// #419): the information adds to the normal matrix, `prior_info (x - x0)` to the gradient, and
/// `L^T` (with `prior_info = L L^T`) to the least-squares rows. Steps are judged on chi-square
/// plus the prior term, but the reported `chi2` is over `astrometry` alone, with
/// `ndof` its number of residual rows; the covariance is the posterior one.
#[allow(clippy::too_many_arguments)]
pub fn fit_orbit_with_prior(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    prior_info: &DMatrix<f64>,
    kernel: &SpiceKernel,
    opts: &FitOptions,
    max_iter: usize,
) -> OrbitFit {
    let opts = FitOptions { nongrav_per_arc: false, ..opts.clone() };
    fit(astrometry, epoch, state, nongrav, None, Some(prior_info), kernel, &opts, max_iter)
}

/// [`fit_orbit`] with per-arc non-gravitational parameters (layup's `per_arc=True`): one state,
/// and the parameters selected in `opts.fit_nongrav` fitted separately for the detections
/// before `epoch` (arc A, starting from `nongrav`) and after it (arc B, starting from
/// `nongrav_arc2`). Place `epoch` between the two apparitions.
#[allow(clippy::too_many_arguments)]
pub fn fit_orbit_per_arc(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    nongrav_arc2: &[f64; 3],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    max_iter: usize,
) -> OrbitFit {
    fit(astrometry, epoch, state, nongrav, Some(*nongrav_arc2), None, kernel, opts, max_iter)
}

#[allow(clippy::too_many_arguments)]
fn fit(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    nongrav_arc2: Option<[f64; 3]>,
    prior_info: Option<&DMatrix<f64>>,
    kernel: &SpiceKernel,
    opts: &FitOptions,
    max_iter: usize,
) -> OrbitFit {
    let active: Vec<usize> = (0..3).filter(|&k| opts.fit_nongrav[k]).collect();
    let per_arc = nongrav_arc2.is_some();
    let nact = active.len();
    let npar = 6 + if per_arc { 2 * nact } else { nact };
    let n = astrometry.len();

    let mut x = *state;
    let mut a = *nongrav;
    let mut b = nongrav_arc2.unwrap_or([0.0; 3]);
    for &k in &active {
        if a[k] == 0.0 {
            a[k] = NONGRAV_SEED;
        }
        if per_arc && b[k] == 0.0 {
            b[k] = NONGRAV_SEED;
        }
    }

    // Prior: its mean (the starting parameters) and the Cholesky factor of its information.
    let pack = |x: &[f64; 6], a: &[f64; 3]| -> DVector<f64> {
        DVector::from_iterator(6 + nact, x.iter().copied().chain(active.iter().map(|&k| a[k])))
    };
    let x0 = pack(&x, &a);
    let prior = match prior_info {
        Some(info) if info.nrows() == npar && info.ncols() == npar => match info.clone().cholesky() {
            Some(ch) => Some((info, ch.l())),
            None => return OrbitFit { epoch, state: *state, ..OrbitFit::failed(FitFlag::PriorNotPositiveDefinite) },
        },
        Some(_) => return OrbitFit { epoch, state: *state, ..OrbitFit::failed(FitFlag::PriorNotPositiveDefinite) },
        None => None,
    };
    let naug = npar + if prior.is_some() { npar } else { 0 };

    // sqrt of the (diagonal) weights, per residual row.
    let w_sqrt: Vec<f64> = (0..n).flat_map(|i| astrometry.rows(i).map(move |k| 1.0 / astrometry.sigma(i, k))).collect();
    let nrows = w_sqrt.len();

    let mut lambda = ARCSEC_PER_RAD * ARCSEC_PER_RAD / 1000.0;
    let mut chi2_prev = f64::INFINITY;
    let mut chi2_final = f64::INFINITY;
    let mut flag = FitFlag::NotConverged;
    let mut normal = DMatrix::<f64>::zeros(npar, npar);
    let mut last_resid: Vec<[f64; 6]> = Vec::new();
    let mut iters = 0;
    let mut accepted_any = false;

    while iters < max_iter {
        let evaluated = if per_arc { residuals_per_arc(astrometry, epoch, &x, &a, &b, kernel, opts, true) } else { residuals(astrometry, epoch, &x, &a, kernel, opts, true) };
        let res = match evaluated {
            Ok(r) => r,
            Err(_) => {
                chi2_final = f64::NAN;
                break;
            }
        };

        // Weighted Jacobian and residuals.
        let mut aw = DMatrix::<f64>::zeros(nrows + naug, npar);
        let mut rhs = DVector::<f64>::zeros(nrows + naug);
        let mut chi2 = 0.0;
        for r in 0..nrows {
            let w = w_sqrt[r];
            for c in 0..npar {
                aw[(r, c)] = w * res.jacobian[r * npar + c];
            }
            let wr = w * res.resid[r];
            rhs[r] = -wr;
            chi2 += wr * wr;
        }
        // Undamped normal matrix B^T W B and gradient B^T W r.
        let bw = aw.rows(0, nrows);
        normal = bw.transpose() * bw;
        let mut grad = -(bw.transpose() * rhs.rows(0, nrows));
        let sl = lambda.sqrt();
        for c in 0..npar {
            aw[(nrows + c, c)] = sl;
        }
        // The prior's rows L^T with right-hand side -L^T (x - x0), and its terms in the normal
        // matrix, the gradient and the objective.
        let mut objective = chi2;
        if let Some((info, l)) = &prior {
            let dx0 = pack(&x, &a) - &x0;
            normal += *info;
            grad += *info * &dx0;
            objective += dx0.dot(&(*info * &dx0));
            let lt = l.transpose();
            let r = -(&lt * &dx0);
            for i in 0..npar {
                for c in 0..npar {
                    aw[(nrows + npar + i, c)] = lt[(i, c)];
                }
                rhs[nrows + npar + i] = r[i];
            }
        }
        let dx = solve_least_squares(aw, rhs);
        last_resid = res.per_detection(n);

        chi2_final = chi2;
        if chi2.is_nan() {
            break;
        }

        let pred = dx.dot(&(lambda * &dx - &grad)).abs();
        let rho = (chi2_prev - objective) / pred;
        if rho > 0.1 {
            lambda *= 0.5;
            for c in 0..6 {
                x[c] += dx[c];
            }
            for (j, &k) in active.iter().enumerate() {
                a[k] += dx[6 + j];
                if per_arc {
                    b[k] += dx[6 + nact + j];
                }
            }
            chi2_prev = objective;
            accepted_any = true;
        } else {
            lambda *= 2.0;
        }

        // Per-parameter tolerance: the absolute `tolerance`, or with `conv_frac` a fraction of the
        // formal sigma from the normal matrix of this iteration (never tighter than `tolerance`).
        let sigma: Vec<f64> = if opts.conv_frac > 0.0 && accepted_any {
            match normal.clone().try_inverse() {
                Some(cov) => (0..npar).map(|i| cov[(i, i)].max(0.0).sqrt()).collect(),
                None => Vec::new(),
            }
        } else {
            Vec::new()
        };
        let tol = |i: usize| match sigma.get(i) {
            Some(&s) if s.is_finite() && s > 0.0 => opts.tolerance.max(opts.conv_frac * s),
            _ => opts.tolerance,
        };
        if dx.iter().enumerate().all(|(i, v)| !v.is_nan() && v.abs() <= tol(i)) {
            flag = FitFlag::Converged;
            break;
        }
        iters += 1;
    }

    // With a prior, the parameters are constrained by it, so the chi-square of these
    // detections is judged over all their rows (layup's convention).
    let ndof = if prior.is_some() { nrows as i64 } else { nrows as i64 - npar as i64 };
    if ndof > 0 && chi2_final / ndof as f64 > opts.chi2_threshold {
        flag = FitFlag::Chi2TooLarge;
    }
    let covariance = normal.clone().try_inverse().unwrap_or_else(|| DMatrix::from_element(npar, npar, f64::NAN));
    if flag == FitFlag::Converged && !active.is_empty() && (6..npar).any(|j| !(covariance[(j, j)] > 0.0) || !covariance[(j, j)].is_finite()) {
        flag = FitFlag::DegenerateCovariance;
    }

    OrbitFit {
        epoch,
        state: x,
        nongrav: a,
        fit_nongrav: opts.fit_nongrav,
        nongrav_arc2: if per_arc { b } else { [0.0; 3] },
        per_arc,
        covariance: covariance.transpose().as_slice().to_vec(),
        npar,
        chi2: chi2_final,
        ndof,
        niter: iters,
        flag,
        residuals: last_resid,
        used: vec![true; astrometry.len()],
    }
}

/// Least-squares solution of `a x = b` (a tall, full column rank) by Householder QR.
fn solve_least_squares(a: DMatrix<f64>, mut b: DVector<f64>) -> DVector<f64> {
    let n = a.ncols();
    let qr = a.qr();
    qr.q_tr_mul(&mut b);
    let r = qr.r();
    let rhs = b.rows(0, n).into_owned();
    r.solve_upper_triangular(&rhs).unwrap_or_else(|| DVector::from_element(n, f64::NAN))
}
