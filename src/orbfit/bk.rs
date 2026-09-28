//! Bernstein–Khushalani linear initial orbit determination (layup's `run_bk_iod`).
//!
//! The orbit is parameterized, following Bernstein & Khushalani (2000), in a tangent-plane frame
//! `(a, b, n0)` around the mean line of sight `n0`: `(alpha, beta)` are the gnomonic coordinates
//! of the direction to the object, `gamma = 1/|r|` (barycentric), and `(adot, bdot, gdot)` are
//! the velocity components along `(a, b, n0)` times `gamma`. Over a short arc and ignoring
//! gravity, each detection's gnomonic coordinates are linear in `(alpha, beta, gamma, adot, bdot)`
//! once `gdot` is fixed at zero:
//!
//! ```text
//! x_i = alpha + adot t_i - gamma (X_i - x_i z_i),   y_i = beta + bdot t_i - gamma (Y_i - y_i z_i)
//! ```
//!
//! with `(X_i, Y_i, z_i)` the observer's position in the frame and `t_i` the light-time corrected
//! time from the epoch. This gives a seed when Gauss's method fails, which it tends to do on short
//! arcs of distant objects, where the three-point geometry is ill-conditioned.

use nalgebra::{Matrix5, Vector3, Vector5};

use crate::data::SPEED_OF_LIGHT;
use crate::orbfit::Astrometry;

/// Right-handed orthonormal frame: `n0` along the mean line of sight, `(a, b)` its tangent plane
/// (`a` from Gram–Schmidt against the z axis, or the x axis when `n0` is near the pole).
fn fiducial(astrometry: &Astrometry) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let mut mean = Vector3::zeros();
    for i in 0..astrometry.len() {
        mean += astrometry.rho_hat(i);
    }
    if mean.norm() < 1e-12 {
        mean = Vector3::x();
    }
    let n0 = mean.normalize();
    let seed = if n0.z.abs() < 0.9 { Vector3::z() } else { Vector3::x() };
    let a = (seed - seed.dot(&n0) * n0).normalize();
    let b = n0.cross(&a);
    (n0, a, b)
}

/// Barycentric J2000 state (AU, AU/day) at TDB Julian date `epoch` from the linear
/// Bernstein–Khushalani fit to `astrometry`, or `None` with fewer than three detections or an
/// unphysical solution (`gamma <= 0` or non-finite).
pub fn bk_iod(astrometry: &Astrometry, epoch: f64) -> Option<[f64; 6]> {
    let n = astrometry.len();
    if n < 3 {
        return None;
    }
    let (n0, a, b) = fiducial(astrometry);

    // Design matrix rows (2 per detection), data and weights.
    let mut rows: Vec<([f64; 5], f64, f64)> = Vec::with_capacity(2 * n);
    let mut ze = Vec::with_capacity(n);
    for i in 0..n {
        let r_obs = astrometry.observer(i);
        let rho = astrometry.rho_hat(i);
        let z = r_obs.dot(&n0);
        let t = (astrometry.epoch[i] + z / SPEED_OF_LIGHT) - epoch;
        let (xo, yo) = (r_obs.dot(&a), r_obs.dot(&b));
        let rn = rho.dot(&n0);
        let (x, y) = (rho.dot(&a) / rn, rho.dot(&b) / rn);
        let (sx, sy) = (astrometry.sigma_ra[i], astrometry.sigma_dec[i]);
        rows.push(([1.0, 0.0, -(xo - x * z), t, 0.0], x, 1.0 / (sx * sx)));
        rows.push(([0.0, 1.0, -(yo - y * z), 0.0, t], y, 1.0 / (sy * sy)));
        ze.push(z);
    }

    // Weighted normal equations H p = g.
    let solve = |w: &dyn Fn(usize) -> f64| -> Vector5<f64> {
        let mut h = Matrix5::zeros();
        let mut g = Vector5::zeros();
        for (k, (row, obs, _)) in rows.iter().enumerate() {
            let wk = w(k);
            for p in 0..5 {
                g[p] += row[p] * wk * obs;
                for q in 0..5 {
                    h[(p, q)] += row[p] * wk * row[q];
                }
            }
        }
        h.col_piv_qr().solve(&g).unwrap_or_else(|| Vector5::from_element(f64::NAN))
    };
    let mut p = solve(&|k| rows[k].2);

    // The relation above is the exact one multiplied through by d_i = 1 - gamma z_i, which
    // reweights each detection by d_i^2: re-solve once with w / d_i^2.
    let gamma0 = p[2];
    if gamma0.is_finite() && gamma0 > 0.0 && ze.iter().all(|z| (1.0 - gamma0 * z) > 0.0) {
        let p2 = solve(&|k| {
            let d = 1.0 - gamma0 * ze[k / 2];
            rows[k].2 / (d * d)
        });
        if p2.iter().all(|v| v.is_finite()) && p2[2] > 0.0 {
            p = p2;
        }
    }

    // The fitted gamma and dots carry a factor s = sqrt(1 + alpha^2 + beta^2).
    let s = (1.0 + p[0] * p[0] + p[1] * p[1]).sqrt();
    let (alpha, beta, gamma, adot, bdot) = (p[0], p[1], p[2] / s, p[3] / s, p[4] / s);
    if ![alpha, beta, gamma, adot, bdot].iter().all(|v| v.is_finite()) || gamma <= 0.0 {
        return None;
    }

    // BK -> Cartesian (gdot = 0).
    let dir = (n0 + alpha * a + beta * b) / s;
    let r = dir / gamma;
    let v = (adot * a + bdot * b) / gamma;
    Some([r.x, r.y, r.z, v.x, v.y, v.z])
}
