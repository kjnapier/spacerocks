//! Gauss's method of initial orbit determination.

use nalgebra::{Matrix3, SMatrix, Vector3};

use crate::data::SPEED_OF_LIGHT;
use crate::orbfit::Astrometry;
use crate::{Observation, SpaceRock};

/// GM of the Sun plus planets (AU^3/day^2), the central mass Gauss's method uses (layup's
/// `GMtotal`).
pub const MU_TOTAL: f64 = 0.0002963092748799319;

/// Candidate orbits from Gauss's method for three detections.
///
/// `rho` are the unit vectors towards the object, `observer` the observer positions and `t` the
/// epochs (TDB Julian dates); the three need not be in time order. Returns one `(epoch, state)`
/// per real root of the eighth-degree range polynomial with `r > min_distance`, largest first:
/// the barycentric state (AU, AU/day) at the middle detection, at an epoch shifted back by
/// `|r| / c` (as in layup, which uses the barycentric distance for this shift).
pub fn gauss_states(rho: [Vector3<f64>; 3], observer: [Vector3<f64>; 3], t: [f64; 3], mu: f64, min_distance: f64) -> Vec<(f64, [f64; 6])> {
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| t[a].total_cmp(&t[b]));
    let [i1, i2, i3] = order;
    let (rho1, rho2, rho3) = (rho[i1], rho[i2], rho[i3]);
    let (o1, o2, o3) = (observer[i1], observer[i2], observer[i3]);
    let (t1, t2, t3) = (t[i1], t[i2], t[i3]);

    let tau1 = t1 - t2;
    let tau3 = t3 - t2;
    let tau = t3 - t1;

    let p1 = rho2.cross(&rho3);
    let p2 = rho1.cross(&rho3);
    let p3 = rho1.cross(&rho2);
    let d0 = rho1.dot(&p1);
    // Coplanar lines of sight leave the ranges indeterminate.
    if d0.abs() < 1e-12 {
        return Vec::new();
    }

    let d = Matrix3::new(
        o1.dot(&p1), o1.dot(&p2), o1.dot(&p3),
        o2.dot(&p1), o2.dot(&p2), o2.dot(&p3),
        o3.dot(&p1), o3.dot(&p2), o3.dot(&p3),
    );

    let a = (1.0 / d0) * (-d[(0, 1)] * (tau3 / tau) + d[(1, 1)] + d[(2, 1)] * (tau1 / tau));
    let b = (1.0 / (6.0 * d0)) * (d[(0, 1)] * ((tau3 * tau3) - (tau * tau)) * (tau3 / tau) + d[(2, 1)] * ((tau * tau) - (tau1 * tau1)) * (tau1 / tau));
    let e = o2.dot(&rho2);
    let o2sq = o2.dot(&o2);

    let aa = -((a * a) + 2.0 * a * e + o2sq);
    let bb = -2.0 * mu * b * (a + e);
    let cc = -mu.powi(2) * (b * b);

    // Companion matrix of r^8 + aa r^6 + bb r^3 + cc.
    let mut m = SMatrix::<f64, 8, 8>::zeros();
    for i in 0..7 {
        m[(i, i + 1)] = 1.0;
    }
    m[(7, 0)] = -cc;
    m[(7, 3)] = -bb;
    m[(7, 6)] = -aa;
    let Some(schur) = m.try_schur(f64::EPSILON, 0) else { return Vec::new() };
    let mut roots: Vec<f64> = schur
        .complex_eigenvalues()
        .iter()
        .filter(|z| z.im.abs() < 1e-10 && z.re > min_distance)
        .map(|z| z.re)
        .collect();
    roots.sort_by(|a, b| b.total_cmp(a));

    roots
        .into_iter()
        .map(|root| {
            let root3 = root.powi(3);
            let num1 = 6.0 * (d[(2, 0)] * (tau1 / tau3) + d[(1, 0)] * (tau / tau3)) * root3 + mu * d[(2, 0)] * ((tau * tau) - (tau1 * tau1)) * (tau1 / tau3);
            let den1 = 6.0 * root3 + mu * ((tau * tau) - (tau3 * tau3));
            let a1 = (1.0 / d0) * ((num1 / den1) - d[(0, 0)]);
            let a2 = a + (mu * b) / root3;
            let num3 = 6.0 * (d[(0, 2)] * (tau3 / tau1) - d[(1, 2)] * (tau / tau1)) * root3 + mu * d[(0, 2)] * ((tau * tau) - (tau3 * tau3)) * (tau3 / tau1);
            let den3 = 6.0 * root3 + mu * ((tau * tau) - (tau1 * tau1));
            let a3 = (1.0 / d0) * ((num3 / den3) - d[(2, 2)]);

            let r1 = o1 + a1 * rho1;
            let r2 = o2 + a2 * rho2;
            let r3 = o3 + a3 * rho3;

            let f1 = 1.0 - 0.5 * (mu / root3) * (tau1 * tau1);
            let f3 = 1.0 - 0.5 * (mu / root3) * (tau3 * tau3);
            let g1 = tau1 - (1.0 / 6.0) * (mu / root3) * (tau1 * tau1 * tau1);
            let g3 = tau3 - (1.0 / 6.0) * (mu / root3) * (tau3 * tau3 * tau3);
            let v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1);

            let epoch = t2 - r2.norm() / SPEED_OF_LIGHT;
            (epoch, [r2.x, r2.y, r2.z, v2.x, v2.y, v2.z])
        })
        .collect()
}

/// Gauss's method on detections `i`, `j` and `k` of `astrometry` (see [`gauss_states`]).
pub fn gauss_astrometry(astrometry: &Astrometry, i: usize, j: usize, k: usize, min_distance: f64) -> Vec<(f64, [f64; 6])> {
    gauss_states(
        [astrometry.rho_hat(i), astrometry.rho_hat(j), astrometry.rho_hat(k)],
        [astrometry.observer(i), astrometry.observer(j), astrometry.observer(k)],
        [astrometry.epoch[i], astrometry.epoch[j], astrometry.epoch[k]],
        MU_TOTAL,
        min_distance,
    )
}

/// Implements Gauss' method for initial orbit determination from three observations.
///
/// Returns the candidate orbits (barycentric J2000 [`SpaceRock`]s at the middle observation,
/// corrected for light time), largest heliocentric distance first, or `None` if there are no
/// real roots beyond `min_distance` (AU).
pub fn gauss(o1: &Observation, o2: &Observation, o3: &Observation, min_distance: f64) -> Option<Vec<SpaceRock>> {
    let a = Astrometry::from_observations(&[o1.clone(), o2.clone(), o3.clone()]).ok()?;
    let rocks: Vec<SpaceRock> = gauss_astrometry(&a, 0, 1, 2, min_distance)
        .into_iter()
        .filter_map(|(epoch, s)| {
            let t = crate::time::Time::new(epoch, "tdb", "jd").ok()?;
            SpaceRock::from_xyz("rock", s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").ok()
        })
        .collect();
    if rocks.is_empty() { None } else { Some(rocks) }
}

/// The triplet of `idx` (indices into `epochs`, in time order) whose outer span is closest to
/// `target_days`, with the middle detection as central as possible while keeping each
/// sub-interval at least 10% of the span. `None` when `idx` already spans less than the target
/// or has fewer than three detections (layup issue #509).
pub fn select_triplet(epochs: &[f64], idx: &[usize], target_days: f64) -> Option<(usize, usize, usize)> {
    let n = idx.len();
    if n < 3 {
        return None;
    }
    let t: Vec<f64> = idx.iter().map(|&i| epochs[i]).collect();
    if t[n - 1] - t[0] <= target_days {
        return None;
    }
    // First index in [lo, hi) with t >= x.
    let bisect_left = |x: f64, lo: usize, hi: usize| lo + t[lo..hi].partition_point(|&v| v < x);
    let mut best = None;
    let mut cost = f64::INFINITY;
    for a in 0..n - 2 {
        let lo = bisect_left(t[a] + target_days, a + 2, n) as isize;
        for c in [lo - 1, lo, lo + 1] {
            if c <= a as isize + 1 || c >= n as isize {
                continue;
            }
            let c = c as usize;
            let span = t[c] - t[a];
            if span <= 0.0 {
                continue;
            }
            let this = (span - target_days).abs();
            if this >= cost {
                continue;
            }
            let floor = 0.1 * span;
            let mid = 0.5 * (t[a] + t[c]);
            let j = bisect_left(mid, a + 1, c) as isize;
            // Candidates ordered by distance from j, lower index first on ties.
            let mut pick = None;
            let mut k = 0isize;
            'outer: loop {
                let mut any = false;
                for b in if k == 0 { vec![j] } else { vec![j - k, j + k] } {
                    if b <= a as isize || b >= c as isize {
                        continue;
                    }
                    any = true;
                    let b = b as usize;
                    if (t[b] - t[a]).min(t[c] - t[b]) >= floor {
                        pick = Some(b);
                        break 'outer;
                    }
                }
                if !any && (j - k <= a as isize && j + k >= c as isize) {
                    break;
                }
                k += 1;
            }
            if let Some(b) = pick {
                best = Some((idx[a], idx[b], idx[c]));
                cost = this;
            }
        }
    }
    best
}

/// Target outer span (days) for [`select_triplet`]: 15 degrees of mean anomaly at an assumed
/// semimajor axis of 2.5 AU.
pub fn gauss_target_days() -> f64 {
    365.25 * 2.5f64.powf(1.5) * 15.0 / 360.0
}
