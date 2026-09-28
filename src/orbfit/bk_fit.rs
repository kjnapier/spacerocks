//! The Bernstein–Khushalani fitting engine: layup's `run_bk_native_fit` (`bk_fit.cpp`, with the
//! parameterization in `bk_basis.cpp`).
//!
//! The orbit is fitted in Bernstein–Khushalani parameters instead of Cartesian ones: the
//! gnomonic direction (alpha, beta) of the barycentric position about a fiducial direction n0,
//! the inverse distance gamma = 1/|r|, and the velocity components along the orthonormal
//! fiducial basis scaled by gamma, (adot, bdot, gdot) = gamma (v·a, v·b, v·n0). The residuals
//! and their Cartesian partials come from the usual variational integration and are carried to
//! these parameters by the chain rule. The line-of-sight velocity gdot, which angles constrain
//! least, gets a fixed Gaussian prior from the bound-orbit condition |v|² < 2 mu / |r|, fixed at
//! the seed: variance 2 mu gamma³ − adot² − bdot² (none if that is not positive). This keeps
//! short arcs of distant objects from running off along the unconstrained direction.

use nalgebra::{DMatrix, DVector, Matrix6, Vector3, Vector6};

use crate::orbfit::residuals::residuals;
use crate::orbfit::{Astrometry, FitFlag, FitOptions, OrbitFit};
use crate::spice::SpiceKernel;

/// layup's MU_SUN (AU³/day²), the mu of the bound-orbit prior.
pub const MU_SUN: f64 = 0.00029591220828559104;

const ARCSEC_PER_RAD: f64 = 206265.0;

/// The fiducial frame: n0 along the mean line of sight, (a, b) spanning its tangent plane.
#[derive(Debug, Clone, Copy)]
pub struct Fiducial {
    pub n0: Vector3<f64>,
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
}

/// layup's `choose_fiducial`, over the detections with optical astrometry.
pub fn choose_fiducial(astrometry: &Astrometry) -> Fiducial {
    let mut mean = Vector3::zeros();
    for i in 0..astrometry.len() {
        if astrometry.has_optical(i) {
            mean += astrometry.rho_hat(i);
        }
    }
    if mean.norm() < 1e-12 {
        mean = Vector3::x();
    }
    let n0 = mean.normalize();
    let seed = if n0.z.abs() < 0.9 { Vector3::z() } else { Vector3::x() };
    let a = (seed - seed.dot(&n0) * n0).normalize();
    let b = n0.cross(&a);
    Fiducial { n0, a, b }
}

/// Bernstein–Khushalani parameters [alpha, beta, gamma, adot, bdot, gdot] of a barycentric
/// Cartesian state (layup's `cartesian_to_bk`).
pub fn cartesian_to_bk(cart: &[f64; 6], f: &Fiducial) -> [f64; 6] {
    let r = Vector3::new(cart[0], cart[1], cart[2]);
    let v = Vector3::new(cart[3], cart[4], cart[5]);
    let gamma = 1.0 / r.norm();
    let rho = gamma * r;
    let u = rho.dot(&f.n0);
    [rho.dot(&f.a) / u, rho.dot(&f.b) / u, gamma, gamma * v.dot(&f.a), gamma * v.dot(&f.b), gamma * v.dot(&f.n0)]
}

/// The line of sight at (alpha, beta) and its derivatives with respect to them.
fn rho_frame(alpha: f64, beta: f64, f: &Fiducial) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let s = (1.0 + alpha * alpha + beta * beta).sqrt();
    let rho = (f.n0 + alpha * f.a + beta * f.b) / s;
    let ra = (f.a - rho.dot(&f.a) * rho) / s;
    let rb = (f.b - rho.dot(&f.b) * rho) / s;
    (rho, ra, rb)
}

/// Barycentric Cartesian state of Bernstein–Khushalani parameters (layup's `bk_to_cartesian`).
pub fn bk_to_cartesian(bk: &[f64; 6], f: &Fiducial) -> [f64; 6] {
    let (rho, _, _) = rho_frame(bk[0], bk[1], f);
    let ig = 1.0 / bk[2];
    let r = ig * rho;
    let v = ig * (bk[3] * f.a + bk[4] * f.b + bk[5] * f.n0);
    [r.x, r.y, r.z, v.x, v.y, v.z]
}

/// d(r, v) / d(alpha, beta, gamma, adot, bdot, gdot) (layup's `dcart_dbk`).
pub fn dcart_dbk(bk: &[f64; 6], f: &Fiducial) -> Matrix6<f64> {
    let (rho, ra, rb) = rho_frame(bk[0], bk[1], f);
    let ig = 1.0 / bk[2];
    let ig2 = ig * ig;
    let mut j = Matrix6::zeros();
    let set = |j: &mut Matrix6<f64>, row: usize, col: usize, v: Vector3<f64>| {
        for k in 0..3 {
            j[(row + k, col)] = v[k];
        }
    };
    set(&mut j, 0, 0, ig * ra);
    set(&mut j, 0, 1, ig * rb);
    set(&mut j, 0, 2, -ig2 * rho);
    set(&mut j, 3, 2, -ig2 * (bk[3] * f.a + bk[4] * f.b + bk[5] * f.n0));
    set(&mut j, 3, 3, ig * f.a);
    set(&mut j, 3, 4, ig * f.b);
    set(&mut j, 3, 5, ig * f.n0);
    j
}

/// Variance of the bound-orbit prior on gdot (layup's `sigma_gdot_sq`): 2 mu gamma³ − adot² −
/// bdot², or infinity (no prior) if that is not positive.
pub fn sigma_gdot_sq(bk: &[f64; 6], mu: f64) -> f64 {
    let rhs = 2.0 * mu * bk[2].powi(3) - bk[3] * bk[3] - bk[4] * bk[4];
    if rhs > 0.0 { rhs } else { f64::INFINITY }
}

/// Fit the gravity-only orbit `state` (barycentric J2000 at TDB `epoch`) to `astrometry` in
/// Bernstein–Khushalani parameters, with the bound-orbit prior on the line-of-sight velocity:
/// layup's `run_bk_native_fit` (`engine="bk_native"`).
///
/// As in layup: the step solves the damped normal equations (`B^T W B + lambda I + P`, with the
/// Jacobian chained to BK parameters) by column-pivoted QR, is accepted when the gain ratio
/// exceeds 0.1, and the fit converges when every component of the BK step is below
/// `opts.tolerance`. The reported chi-square includes the prior's term, ndof is rows − 6, and the
/// covariance is carried back to Cartesian. Unlike layup (which uses only RA and Dec), rate and
/// radar rows are fitted too, through the same chain rule.
pub fn fit_orbit_bk(astrometry: &Astrometry, epoch: f64, state: &[f64; 6], kernel: &SpiceKernel, opts: &FitOptions, max_iter: usize) -> OrbitFit {
    let n = astrometry.len();
    let mut out = OrbitFit { epoch, state: *state, covariance: vec![0.0; 36], ..OrbitFit::failed(FitFlag::NotConverged) };
    if n < 3 {
        return out;
    }
    let opts = FitOptions { fit_nongrav: [false; 3], nongrav_per_arc: false, ..opts.clone() };
    let fid = choose_fiducial(astrometry);
    let mut bk = cartesian_to_bk(state, &fid);
    let sgsq = sigma_gdot_sq(&bk, MU_SUN);
    let mut prior = Matrix6::<f64>::zeros();
    prior[(5, 5)] = if sgsq.is_finite() && sgsq > 0.0 { 1.0 / sgsq } else { 0.0 };

    let weights: Vec<f64> = (0..n).flat_map(|i| astrometry.rows(i).map(move |k| astrometry.sigma(i, k).powi(-2))).collect();
    let nrows = weights.len();

    let mut lambda = ARCSEC_PER_RAD * ARCSEC_PER_RAD / 1000.0;
    let (mut chi2_prev, mut chi2_cur) = (f64::INFINITY, f64::INFINITY);
    let mut c_final: Option<Matrix6<f64>> = None;
    let mut last_resid: Vec<[f64; 6]> = Vec::new();
    let mut flag = FitFlag::NotConverged;
    let mut iters = 0;
    while iters < max_iter {
        let cart = bk_to_cartesian(&bk, &fid);
        let res = match residuals(astrometry, epoch, &cart, &[0.0; 3], kernel, &opts, true) {
            Ok(r) => r,
            Err(_) => {
                chi2_cur = f64::NAN;
                break;
            }
        };
        last_resid = res.per_detection(n);
        let b_cart = DMatrix::from_row_slice(nrows, 6, &res.jacobian);
        let j = dcart_dbk(&bk, &fid);
        let b_bk = &b_cart * DMatrix::from_column_slice(6, 6, j.as_slice());
        let r = DVector::from_column_slice(&res.resid);
        let w = DVector::from_column_slice(&weights);
        let btw = DMatrix::from_fn(6, nrows, |p, q| b_bk[(q, p)] * w[q]);
        let c_data = Matrix6::from_iterator((&btw * &b_bk).iter().copied());
        let bk_vec = Vector6::from_column_slice(&bk);
        let grad = Vector6::from_iterator((&btw * &r).iter().copied()) + prior * bk_vec;
        let chi2_data: f64 = (0..nrows).map(|q| r[q] * r[q] * w[q]).sum();
        chi2_cur = chi2_data + bk_vec.dot(&(prior * bk_vec));
        if c_final.is_none() {
            c_final = Some(c_data + prior); // the seed's, in case no step is accepted
        }

        let c = c_data + Matrix6::identity() * lambda + prior;
        let dx = DMatrix::from_column_slice(6, 6, c.as_slice()).col_piv_qr().solve(&DVector::from_column_slice((-grad).as_slice()));
        let dx = match dx {
            Some(d) => Vector6::from_column_slice(d.as_slice()),
            None => Vector6::from_element(f64::NAN),
        };
        let den = dx.dot(&(lambda * dx - grad)).abs();
        let rho = if den > 0.0 { (chi2_prev - chi2_cur) / den } else { -1.0 };
        if rho > 0.1 {
            lambda *= 0.5;
            for k in 0..6 {
                bk[k] += dx[k];
            }
            chi2_prev = chi2_cur;
            c_final = Some(c_data + prior);
        } else {
            lambda *= 2.0;
        }
        if !chi2_cur.is_nan() && dx.iter().all(|v| !v.is_nan() && v.abs() <= opts.tolerance) {
            flag = FitFlag::Converged;
            break;
        }
        iters += 1;
    }

    let ndof = nrows as i64 - 6;
    if flag == FitFlag::Converged && ndof > 0 && chi2_cur / ndof as f64 > opts.chi2_threshold {
        flag = FitFlag::Chi2TooLarge;
    }
    if let Some(c) = c_final.and_then(|c| c.try_inverse()) {
        let j = dcart_dbk(&bk, &fid);
        let cov = j * c * j.transpose();
        out.covariance = (0..36).map(|q| cov[(q / 6, q % 6)]).collect();
    }
    out.state = bk_to_cartesian(&bk, &fid);
    out.chi2 = chi2_cur;
    out.ndof = ndof;
    out.niter = iters;
    out.flag = flag;
    out.residuals = last_resid;
    out.used = vec![true; n];
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fid() -> Fiducial {
        let n0 = Vector3::new(0.3, -0.8, 0.52).normalize();
        let z: Vector3<f64> = Vector3::z();
        let a = (z - z.dot(&n0) * n0).normalize();
        Fiducial { n0, a, b: n0.cross(&a) }
    }

    #[test]
    fn round_trip_and_jacobian() {
        let f = fid();
        let cart = [12.0, -30.0, 18.0, 0.0021, 0.0009, -0.0012];
        let bk = cartesian_to_bk(&cart, &f);
        let back = bk_to_cartesian(&bk, &f);
        for k in 0..6 {
            assert!((back[k] - cart[k]).abs() < 1e-12 * cart[k].abs().max(1e-3), "{} {} {}", k, back[k], cart[k]);
        }
        let j = dcart_dbk(&bk, &f);
        for c in 0..6 {
            let h = 1e-7 * bk[c].abs().max(1e-6);
            let (mut p, mut m) = (bk, bk);
            p[c] += h;
            m[c] -= h;
            let (cp, cm) = (bk_to_cartesian(&p, &f), bk_to_cartesian(&m, &f));
            for r in 0..6 {
                let fd = (cp[r] - cm[r]) / (2.0 * h);
                assert!((fd - j[(r, c)]).abs() <= 1e-6 * j.column(c).amax().max(1e-12), "J[{},{}]: {} vs {}", r, c, fd, j[(r, c)]);
            }
        }
    }
}
