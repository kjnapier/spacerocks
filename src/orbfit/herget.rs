//! Herget's method of initial orbit determination, as layup implements it
//! (`utilities/herget_iod.py`, with its universal-variable Kepler propagator
//! `utilities/universal_kepler.py`).
//!
//! The ranges to the first and last detections of an arc, rho_1 and rho_n, are adjusted until the
//! orbit through the two positions (found by shooting with two-body motion, then followed with
//! ASSIST's force model) fits the detections in between. Each iteration linearizes the residuals in
//! the two ranges and solves the 2x2 normal equations.

use nalgebra::{Matrix3, Vector3};

use crate::orbfit::residuals::trial_simulation;
use crate::orbfit::{Astrometry, FitOptions};
use crate::spice::SpiceKernel;

/// layup's MU_SUN and SPEED_OF_LIGHT (AU^3/day^2, AU/day).
const MU_SUN: f64 = 0.00029591220828559104;
const SPEED_OF_LIGHT: f64 = 173.1446326742403;

// ---- Universal-variable Kepler propagation (Danby 1988), layup's port of universal-kepler.c ----

const NEWTON_MAX: usize = 6;
const LAGCON_MAX: usize = 15;
const TOL: f64 = 1e-13;
const F_FLOOR: f64 = 1e-8;
const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

fn converged(ds: f64, s: f64) -> bool {
    ds.abs() <= TOL * 1f64.max(s.abs())
}

/// Stumpff functions c0..c5 by Mikkola's argument four-folding.
fn stumpff_c(z: f64) -> [f64; 6] {
    let mut h = z;
    let mut k = 0;
    while h.abs() >= 0.1 {
        h *= 0.25;
        k += 1;
    }
    let mut c4 = (1.0 - h * (1.0 - h * (1.0 - h / 90.0 / (1.0 + h / 132.0)) / 56.0) / 30.0) / 24.0;
    let mut c5 = (1.0 - h * (1.0 - h * (1.0 - h / 110.0 / (1.0 + h / 156.0)) / 72.0) / 42.0) / 120.0;
    for _ in 0..k {
        let c3 = 1.0 / 6.0 - h * c5;
        let c2 = 0.5 - h * c4;
        c5 = (c5 + c4 + c2 * c3) / 16.0;
        c4 = c3 * (2.0 - h * c3) / 8.0;
        h *= 4.0;
    }
    let c3 = 1.0 / 6.0 - z * c5;
    let c2 = 0.5 - z * c4;
    let c1 = 1.0 - z * c3;
    let c0 = 1.0 - z * c2;
    [c0, c1, c2, c3, c4, c5]
}

fn initial_guess(gm: f64, dt: f64, r0: f64, alpha: f64, u: f64) -> f64 {
    if (dt / r0).abs() <= 0.2 {
        return dt / r0 - (dt * dt * u) / (2.0 * r0.powi(3));
    }
    if alpha <= 0.0 {
        let a = gm / alpha;
        let en = (-gm / (a * a * a)).sqrt();
        let ch = 1.0 - r0 / a;
        let sh = u / (-a * gm).sqrt();
        let e = (ch * ch - sh * sh).sqrt();
        let dm = en * dt;
        return if dm < 0.0 {
            -((-2.0 * dm + 1.8 * e) / (ch - sh)).ln() / (-alpha).sqrt()
        } else {
            ((2.0 * dm + 1.8 * e) / (ch + sh)).ln() / (-alpha).sqrt()
        };
    }
    let a = gm / alpha;
    let en = (gm / (a * a * a)).sqrt();
    let ec = 1.0 - r0 / a;
    let es = u / (en * a * a);
    let n_rev = (en * dt / TWO_PI).trunc();
    let dt_red = dt - n_rev * TWO_PI / en;
    let y = en * dt_red - es;
    let (xx, mut yy) = (ec, es);
    let h = en * dt_red;
    let mut omx = h / (1.0 - xx);
    let (k0x, k0y) = (-yy * omx, xx * omx);
    let (mut xx1, mut yy1) = (xx + k0x / 2.0, yy + k0y / 2.0);
    omx = h / (1.0 - xx1);
    let (k1x, k1y) = (-yy1 * omx, xx1 * omx);
    xx1 = xx + k1x / 2.0;
    yy1 = yy + k1y / 2.0;
    omx = h / (1.0 - xx1);
    let (k2x, k2y) = (-yy1 * omx, xx1 * omx);
    xx1 = xx + k2x;
    omx = h / (1.0 - xx1);
    let k3y = xx1 * omx;
    yy += (k0y + 2.0 * (k1y + k2y) + k3y) / 6.0;
    let root_alpha = alpha.sqrt();
    (y + yy) / root_alpha + n_rev * TWO_PI / root_alpha
}

/// Propagate `state` by `dt` under two-body motion (`gm`), optionally carrying a deviation of
/// the initial state along: layup's `universal_step`.
fn universal_step(gm: f64, dt: f64, state: &[f64; 6], variation: Option<&[f64; 6]>) -> Result<([f64; 6], Option<[f64; 6]>), String> {
    let r0v = Vector3::new(state[0], state[1], state[2]);
    let v0v = Vector3::new(state[3], state[4], state[5]);
    let r0 = r0v.norm();
    if r0 == 0.0 {
        return Err("zero position".into());
    }
    let v0s = v0v.dot(&v0v);
    let u = r0v.dot(&v0v);
    let alpha = 2.0 * gm / r0 - v0s;
    let zeta = gm - alpha * r0;

    let fvals = |s: f64| {
        let c = stumpff_c(s * s * alpha);
        let (c0, c1, c2, c3) = (c[0], c[1] * s, c[2] * s * s, c[3] * s * s * s);
        let f = r0 * c1 + u * c2 + gm * c3 - dt;
        let fp = r0 * c0 + u * c1 + gm * c2;
        let fpp = zeta * c1 + u * c0;
        let fppp = zeta * c0 - u * alpha * c1;
        (f, fp, fpp, fppp)
    };
    let s_guess = initial_guess(gm, dt, r0, alpha, u);
    let mut s = s_guess;
    let mut ds = f64::INFINITY;
    for _ in 0..NEWTON_MAX {
        let (f, fp, fpp, fppp) = fvals(s);
        ds = -f / fp;
        ds = -f / (fp + ds * fpp / 2.0);
        ds = -f / (fp + ds * fpp / 2.0 + ds * ds * fppp / 6.0);
        s += ds;
        if converged(ds, s) {
            break;
        }
    }
    if !converged(ds, s) {
        s = s_guess;
        let ln = 5.0;
        for _ in 0..LAGCON_MAX {
            let (f, fp, fpp, _) = fvals(s);
            let disc = (ln - 1.0) * (ln - 1.0) * fp * fp - (ln - 1.0) * ln * f * fpp;
            ds = -ln * f / (fp + disc.abs().sqrt().copysign(fp));
            s += ds;
            if converged(ds, s) {
                break;
            }
        }
        if !converged(ds, s) {
            return Err(format!("Kepler equation did not converge (dt {}, r0 {}, alpha {})", dt, r0, alpha));
        }
    }
    let c = stumpff_c(s * s * alpha);
    let (g0, g1, g2, g3, g4, g5) = (c[0], c[1] * s, c[2] * s * s, c[3] * s.powi(3), c[4] * s.powi(4), c[5] * s.powi(5));
    let r = r0 * g0 + u * g1 + gm * g2;

    let f = 1.0 - (gm / r0) * g2;
    let g = dt - gm * g3;
    let fdot = -(gm / (r * r0)) * g1;
    let gdot = if f.abs() > F_FLOOR { (1.0 + g * fdot) / f } else { 1.0 - (gm / r) * g2 };
    let pos = r0v * f + v0v * g;
    let vel = r0v * fdot + v0v * gdot;
    let out = [pos.x, pos.y, pos.z, vel.x, vel.y, vel.z];

    let var_out = variation.map(|d| {
        let dr = Vector3::new(d[0], d[1], d[2]);
        let dv = Vector3::new(d[3], d[4], d[5]);
        let r0pr = r0v.dot(&dr) / r0;
        let alphapr = -(2.0 * gm / (r0 * r0)) * r0pr - 2.0 * v0v.dot(&dv);
        let upr = r0v.dot(&dv) + v0v.dot(&dr);
        let zetapr = -alpha * r0pr - r0 * alphapr;
        let g1a = 0.5 * (g3 - s * g2);
        let g2a = 0.5 * (2.0 * g4 - s * g3);
        let g3a = 0.5 * (3.0 * g5 - s * g4);
        let spr = -(s * r0pr + g3 * zetapr + g2 * upr + (g3a * zeta + u * g2a) * alphapr) / r;
        let g1pr = g0 * spr + g1a * alphapr;
        let g2pr = g1 * spr + g2a * alphapr;
        let g3pr = g2 * spr + g3a * alphapr;
        let rpr = r0pr + g1 * upr + g2 * zetapr + u * g1pr + zeta * g2pr;
        let fpr = (gm * g2 / (r0 * r0)) * r0pr - (gm / r0) * g2pr;
        let gpr = -gm * g3pr;
        let fdotpr = (gm / (r * r * r0)) * g1 * rpr + (gm / (r * r0 * r0)) * g1 * r0pr - (gm / (r * r0)) * g1pr;
        let gdotpr = (gm / (r * r)) * g2 * rpr - (gm / r) * g2pr;
        let p = dr * f + dv * g + r0v * fpr + v0v * gpr;
        let v = dr * fdot + dv * gdot + r0v * fdotpr + v0v * gdotpr;
        [p.x, p.y, p.z, v.x, v.y, v.z]
    });
    Ok((out, var_out))
}

/// The 6x6 two-body state transition matrix d state(t0 + dt) / d state(t0).
fn state_transition_matrix(gm: f64, dt: f64, state: &[f64; 6]) -> Result<[[f64; 6]; 6], String> {
    let mut stm = [[0.0; 6]; 6];
    for j in 0..6 {
        let mut e = [0.0; 6];
        e[j] = 1.0;
        let (_, v) = universal_step(gm, dt, state, Some(&e))?;
        let v = v.unwrap();
        for i in 0..6 {
            stm[i][j] = v[i];
        }
    }
    Ok(stm)
}

// ---- Herget -------------------------------------------------------------------------------------

/// Velocity at t1 that carries position r1 to rn at tn (two-body, shooting one velocity component
/// at a time), and the velocity it arrives with: layup's `find_velocity`.
fn find_velocity(t1: f64, tn: f64, r1: &Vector3<f64>, rn: &Vector3<f64>, tolerance: f64, max_iter: usize) -> Result<(Vector3<f64>, Vector3<f64>), String> {
    let dt = tn - t1;
    let mut s1 = [r1.x, r1.y, r1.z, (rn.x - r1.x) / dt, (rn.y - r1.y) / dt, (rn.z - r1.z) / dt];
    let mut sn = [rn.x + tolerance.abs() + 1.0, rn.y + tolerance.abs() + 1.0, rn.z + tolerance.abs() + 1.0, 0.0, 0.0, 0.0];
    let mut iter = 0;
    while (Vector3::new(sn[0], sn[1], sn[2]) - rn).norm() > tolerance && iter < max_iter {
        // layup's find_new_vel_with_universal_kepler
        for i in 0..3 {
            let mut e = [0.0; 6];
            e[i + 3] = 1.0;
            let (st, var) = universal_step(MU_SUN, dt, &s1, Some(&e))?;
            let var = var.unwrap();
            let q = Vector3::new(st[0], st[1], st[2]);
            let rv = q + Vector3::new(var[0], var[1], var[2]);
            let mag = (rv - q).dot(&(rn - q)) / (rv - q).dot(&(rv - q));
            s1[i + 3] += mag;
            sn = st;
        }
        iter += 1;
    }
    Ok((Vector3::new(s1[3], s1[4], s1[5]), Vector3::new(sn[3], sn[4], sn[5])))
}

fn inverse(m: Matrix3<f64>) -> Result<Matrix3<f64>, String> {
    m.try_inverse().ok_or_else(|| "singular state transition block".to_string())
}

fn block(stm: &[[f64; 6]; 6], r0: usize, c0: usize) -> Matrix3<f64> {
    Matrix3::from_fn(|i, j| stm[r0 + i][c0 + j])
}

/// One Herget correction (layup's `find_drho`): returns (delta rho_1, delta rho_n, the state at
/// t1, the normalized determinant of the normal equations).
#[allow(clippy::too_many_arguments)]
fn find_drho(
    a: &Astrometry,
    idx: &[usize],
    epochs: &[f64],
    t1: f64,
    tn: f64,
    r1: &Vector3<f64>,
    rn: &Vector3<f64>,
    tolerance: f64,
    rho_1: f64,
    rho_hat_1: &Vector3<f64>,
    rho_hat_n: &Vector3<f64>,
    max_iter: usize,
    kernel: &SpiceKernel,
    opts: &FitOptions,
) -> Result<(f64, f64, [f64; 6], f64), String> {
    let (v1, vn) = find_velocity(t1, tn, r1, rn, tolerance * rho_1 / 100.0, max_iter)?;
    let s1 = [r1.x, r1.y, r1.z, v1.x, v1.y, v1.z];
    let phi = state_transition_matrix(MU_SUN, tn - t1, &s1)?;
    let inv_rv = inverse(block(&phi, 0, 3))?;
    let var_v1 = -(inv_rv * block(&phi, 0, 0)) * rho_hat_1;
    let var_vn = (block(&phi, 3, 3) * inv_rv) * rho_hat_n;
    let sn = [rn.x, rn.y, rn.z, vn.x, vn.y, vn.z];

    // Predicted unit vectors along the orbit, and along the orbit displaced by the full
    // (unit-range) variation, both projected on each detection's tangent basis.
    let project = |start: f64, state: &[f64; 6], var: &Vector3<f64>, dir: &Vector3<f64>| -> Result<Vec<(f64, f64, f64, f64)>, String> {
        let mut sim = trial_simulation(start, state, &[0.0; 3], &[], kernel, opts, false).map_err(|e| e.to_string())?;
        sim.add_variation_from_state("herget", [dir.x, dir.y, dir.z, var.x, var.y, var.z], "rock").map_err(|e| e.to_string())?;
        let jd_ref = sim.state.jd_ref;
        let mut out = Vec::with_capacity(idx.len());
        for (k, &i) in idx.iter().enumerate() {
            sim.integrate_rel(epochs[k] - jd_ref, kernel).map_err(|e| e.to_string())?;
            let r = sim.state.particles[0].position;
            let dr = sim.state.variational_particles[0].position;
            let rho = r - a.observer(i);
            let (av, dv) = a.tangent_basis(i);
            let u = rho / rho.norm();
            let w = (rho + dr) / (rho + dr).norm();
            out.push((u.dot(&av), u.dot(&dv), w.dot(&av), w.dot(&dv)));
        }
        Ok(out)
    };
    let p1 = project(t1, &s1, &var_v1, rho_hat_1)?;
    let pn = project(tn, &sn, &var_vn, rho_hat_n)?;

    let (mut a1b, mut a2b, mut a1sq, mut a2sq, mut a1a2) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (x, y) in p1.iter().zip(&pn) {
        for (b, w1, w2) in [(x.0, x.2, y.2), (x.1, x.3, y.3)] {
            let (c1, c2) = (b - w1, b - w2);
            a1b += c1 * b;
            a2b += c2 * b;
            a1sq += c1 * c1;
            a2sq += c2 * c2;
            a1a2 += c1 * c2;
        }
    }
    let d1 = (a1b * a2sq - a2b * a1a2) / (a1a2 * a1a2 - a1sq * a2sq);
    let dn = -(a2b + d1 * a1a2) / a2sq;
    let det = (a1sq * a2sq - a1a2 * a1a2) / (a1sq * a2sq);
    Ok((d1, dn, s1, det))
}

/// Herget's method on the detections `idx` of `astrometry` (in time order), starting from the
/// range `initial_rho` (AU) to the first and last: layup's `herget_with_assist` (tolerance 0.003
/// AU, 100 iterations). Returns the epoch of the first detection and the state there, or `None`
/// if the ranges don't settle or the normal equations become degenerate.
pub fn herget(astrometry: &Astrometry, idx: &[usize], initial_rho: f64, kernel: &SpiceKernel, opts: &FitOptions) -> Option<(f64, [f64; 6])> {
    const TOLERANCE: f64 = 0.003;
    const MAX_ITER: usize = 100;
    if idx.len() < 3 {
        return None;
    }
    let (first, last) = (idx[0], idx[idx.len() - 1]);
    let (re1, ren) = (astrometry.observer(first), astrometry.observer(last));
    let (rh1, rhn) = (astrometry.rho_hat(first), astrometry.rho_hat(last));
    let (t1, tn) = (astrometry.epoch[first], astrometry.epoch[last]);
    let (mut rho_1, mut rho_n) = (initial_rho, initial_rho);
    let mut r1 = re1 + rh1 * rho_1;
    let mut rn = ren + rhn * rho_n;
    let (mut d1, mut dn) = (TOLERANCE + 1.0, TOLERANCE + 1.0);
    let mut dets: Vec<f64> = vec![1.0];
    let mut state = [f64::NAN; 6];
    let mut iteration = 0;
    let median = |v: &Vec<f64>| {
        let mut s = v.clone();
        s.sort_by(|a, b| a.total_cmp(b));
        let m = s.len();
        if m % 2 == 1 { s[m / 2] } else { 0.5 * (s[m / 2 - 1] + s[m / 2]) }
    };
    while (d1.abs() + dn.abs()) / 2.0 > TOLERANCE && iteration < MAX_ITER && median(&dets) > 1e-3 {
        let lt = (rho_1 + rho_n) / (2.0 * SPEED_OF_LIGHT);
        let epochs: Vec<f64> = idx.iter().map(|&i| astrometry.epoch[i] - lt).collect();
        let (a, b, s1, det) = find_drho(astrometry, idx, &epochs, t1, tn, &r1, &rn, TOLERANCE, rho_1, &rh1, &rhn, MAX_ITER, kernel, opts).ok()?;
        if !(a.is_finite() && b.is_finite() && det.is_finite()) {
            return None;
        }
        d1 = if a.abs() > rho_1 / 2.0 { a.signum() * rho_1 / 2.0 } else { a };
        dn = if b.abs() > rho_n / 2.0 { b.signum() * rho_n / 2.0 } else { b };
        if iteration == 0 {
            dets = vec![det];
        } else {
            dets.push(det);
        }
        state = s1;
        rho_1 -= d1;
        r1 = re1 + rh1 * rho_1;
        rho_n -= dn;
        rn = ren + rhn * rho_n;
        iteration += 1;
    }
    if iteration >= MAX_ITER || median(&dets) <= 1e-3 || !state.iter().all(|v| v.is_finite()) {
        return None;
    }
    Some((t1, state))
}

/// layup's `herget_iod`: Herget on the detections `idx`, from initial ranges of 2, then 5, then
/// 40 AU, until one works.
pub fn herget_iod(astrometry: &Astrometry, idx: &[usize], kernel: &SpiceKernel, opts: &FitOptions) -> Option<(f64, [f64; 6])> {
    [2.0, 5.0, 40.0].iter().find_map(|&rho| herget(astrometry, idx, rho, kernel, opts))
}
