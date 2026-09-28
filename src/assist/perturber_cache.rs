//! Precomputed perturber ephemeris for IAS15.
//!
//! Every IAS15 step needs the barycentric states of all perturbers at eight sub-step times.
//! Reading them from SPK files means walking each body's segment chain (for example
//! Moon -> Earth-Moon barycenter -> SSB) and evaluating several Chebyshev records per body.
//! [`PerturberCache`] instead fits each body's barycentric state once, over a fixed time span,
//! with one Chebyshev polynomial per interval, and evaluates those fits directly.
//!
//! The fits are built adaptively: each body starts with long intervals, which are halved until
//! the fit matches the kernel to within [`POSITION_TOLERANCE`] (AU) and
//! [`VELOCITY_TOLERANCE`] (AU/day) at test points between the fitting nodes, so the cache
//! reproduces the kernel far below any physically meaningful level. A cache is immutable and
//! can be shared between simulations and threads (`Arc<PerturberCache>`). Epochs outside the
//! cached span fall back to the kernel.

use crate::spice::error::Result;
use crate::spice::{SpiceError, SpiceKernel};

/// Number of Chebyshev coefficients per component per interval (degree 15).
const NCOEF: usize = 16;
/// Longest and shortest interval tried (days).
const MAX_INTERVAL: f64 = 32.0;
const MIN_INTERVAL: f64 = 0.125;
/// Maximum allowed fit error at the test points, on top of the kernel's own evaluation noise
/// (see [`EPOCH_NOISE`]).
pub const POSITION_TOLERANCE: f64 = 1e-14; // AU (~1.5 mm)
pub const VELOCITY_TOLERANCE: f64 = 1e-14; // AU/day (~2e-8 m/s)
/// Epochs reach the kernel as ephemeris seconds past J2000 computed from a Julian date, which
/// rounds them to ~1e-7 s; the kernel's output therefore jitters by about velocity x 1e-12 days
/// (4e-14 AU for Mercury). No fit can do better, so the tolerance allows for that.
pub const EPOCH_NOISE: f64 = 4e-12; // days
/// Test points (fractions of an interval) used to validate a fit, between the nodes and near
/// both ends.
const TEST_POINTS: [f64; 8] = [0.0005, 0.047, 0.19, 0.36, 0.53, 0.71, 0.88, 0.9995];

#[derive(Debug, Clone)]
struct BodyFit {
    interval: f64,
    n_intervals: usize,
    /// `[interval][coefficient][component (x, y, z, vx, vy, vz)]`
    coeffs: Vec<f64>,
}

/// Chebyshev fits of perturber barycentric J2000 states (AU, AU/day) over a time span.
#[derive(Debug, Clone)]
pub struct PerturberCache {
    ids: Vec<i32>,
    t_start: f64,
    t_end: f64,
    fits: Vec<BodyFit>,
}

impl PerturberCache {
    /// Fit the barycentric states of `ids` over `[t_start, t_end]` (TDB Julian dates).
    /// Fails if the kernel does not cover the whole span.
    pub fn build(kernel: &SpiceKernel, ids: &[i32], t_start: f64, t_end: f64) -> Result<PerturberCache> {
        if !(t_end > t_start) {
            return Err(SpiceError::Config(format!("invalid cache span [{}, {}]", t_start, t_end)));
        }
        let span = t_end - t_start;
        let mut fits: Vec<Option<BodyFit>> = vec![None; ids.len()];
        let mut interval = MAX_INTERVAL.min(span);

        // Chebyshev nodes on [-1, 1] and the cosine table for the fit.
        let nodes: Vec<f64> = (0..NCOEF)
            .map(|k| (std::f64::consts::PI * (k as f64 + 0.5) / NCOEF as f64).cos())
            .collect();

        loop {
            let pending: Vec<usize> = (0..ids.len()).filter(|&b| fits[b].is_none()).collect();
            if pending.is_empty() {
                break;
            }
            let pending_ids: Vec<i32> = pending.iter().map(|&b| ids[b]).collect();
            let n_intervals = (span / interval).ceil().max(1.0) as usize;
            let interval_len = span / n_intervals as f64;

            // Sample every pending body at the nodes of every interval.
            let mut samples = vec![[0.0f64; 6]; pending.len()];
            let mut coeffs: Vec<Vec<f64>> = vec![vec![0.0; n_intervals * 6 * NCOEF]; pending.len()];
            let mut values = vec![vec![[0.0f64; 6]; NCOEF]; pending.len()];
            for k in 0..n_intervals {
                let a = t_start + k as f64 * interval_len;
                // Epochs are Julian dates, whose resolution (~40 us) is coarse compared with
                // the fit accuracy, so the sampled epochs are not exactly at the Chebyshev
                // nodes. Interpolate at the local coordinates of the epochs actually sampled,
                // computed exactly as `cheb_eval` computes them.
                let mut vander = nalgebra::DMatrix::<f64>::zeros(NCOEF, NCOEF);
                for (n, &x) in nodes.iter().enumerate() {
                    let t = a + 0.5 * (x + 1.0) * interval_len;
                    kernel.barycentric_states_au(&pending_ids, t, &mut samples)?;
                    for p in 0..pending.len() {
                        values[p][n] = samples[p];
                    }
                    let xa = local_coordinate(t, t_start, interval_len, n_intervals).1;
                    let (mut t0, mut t1) = (1.0, xa);
                    vander[(n, 0)] = 1.0;
                    vander[(n, 1)] = xa;
                    for j in 2..NCOEF {
                        let t2 = 2.0 * xa * t1 - t0;
                        vander[(n, j)] = t2;
                        t0 = t1;
                        t1 = t2;
                    }
                }
                let lu = vander.lu();
                let mut rhs = nalgebra::DMatrix::<f64>::zeros(NCOEF, pending.len() * 6);
                for p in 0..pending.len() {
                    for c in 0..6 {
                        for n in 0..NCOEF {
                            rhs[(n, p * 6 + c)] = values[p][n][c];
                        }
                    }
                }
                let sol = lu.solve(&rhs).ok_or_else(|| SpiceError::Config("singular Chebyshev fit".into()))?;
                for p in 0..pending.len() {
                    for c in 0..6 {
                        for j in 0..NCOEF {
                            coeffs[p][k * 6 * NCOEF + j * 6 + c] = sol[(j, p * 6 + c)];
                        }
                    }
                }
            }

            // Validate between the nodes.
            let mut max_err = vec![(0.0f64, 0.0f64); pending.len()];
            // Largest speed and (centripetal) acceleration of each body over the span.
            let mut motion_scale = vec![(0.0f64, 0.0f64); pending.len()];
            let mut fit_out = [0.0f64; 6];
            for k in 0..n_intervals {
                let a = t_start + k as f64 * interval_len;
                for &f in TEST_POINTS.iter() {
                    let t = a + f * interval_len;
                    kernel.barycentric_states_au(&pending_ids, t, &mut samples)?;
                    for p in 0..pending.len() {
                        cheb_eval(&coeffs[p], t_start, interval_len, n_intervals, t, &mut fit_out);
                        let dp = ((fit_out[0] - samples[p][0]).powi(2)
                            + (fit_out[1] - samples[p][1]).powi(2)
                            + (fit_out[2] - samples[p][2]).powi(2))
                        .sqrt();
                        let dv = ((fit_out[3] - samples[p][3]).powi(2)
                            + (fit_out[4] - samples[p][4]).powi(2)
                            + (fit_out[5] - samples[p][5]).powi(2))
                        .sqrt();
                        let r = (samples[p][0].powi(2) + samples[p][1].powi(2) + samples[p][2].powi(2)).sqrt();
                        let v2 = samples[p][3].powi(2) + samples[p][4].powi(2) + samples[p][5].powi(2);
                        motion_scale[p].0 = motion_scale[p].0.max(v2.sqrt());
                        motion_scale[p].1 = motion_scale[p].1.max(v2 / r.max(1e-3));
                        max_err[p].0 = max_err[p].0.max(dp);
                        max_err[p].1 = max_err[p].1.max(dv);
                    }
                }
            }

            let last_level = interval <= MIN_INTERVAL;
            for (p, &b) in pending.iter().enumerate() {
                let (speed, accel) = motion_scale[p];
                let ok = max_err[p].0 <= POSITION_TOLERANCE + EPOCH_NOISE * speed
                    && max_err[p].1 <= VELOCITY_TOLERANCE + EPOCH_NOISE * accel;
                if ok {
                    fits[b] = Some(BodyFit {
                        interval: interval_len,
                        n_intervals,
                        coeffs: std::mem::take(&mut coeffs[p]),
                    });
                } else if last_level {
                    return Err(SpiceError::Config(format!(
                        "could not fit body {} to tolerance (position error {:e} AU, velocity error {:e} AU/day)",
                        ids[b], max_err[p].0, max_err[p].1
                    )));
                }
            }
            interval /= 2.0;
        }

        Ok(PerturberCache {
            ids: ids.to_vec(),
            t_start,
            t_end,
            fits: fits.into_iter().map(|f| f.unwrap()).collect(),
        })
    }

    /// NAIF IDs of the cached bodies, in order.
    pub fn ids(&self) -> &[i32] {
        &self.ids
    }

    /// Cached span (TDB Julian dates).
    pub fn span(&self) -> (f64, f64) {
        (self.t_start, self.t_end)
    }

    /// Interval length (days) used for each body.
    pub fn intervals(&self) -> Vec<(i32, f64)> {
        self.ids.iter().zip(&self.fits).map(|(&id, f)| (id, f.interval)).collect()
    }

    /// Fill `out` with the states of `ids` at `t`. Returns false (and leaves `out` unspecified)
    /// if `ids` is not the cached body list or `t` is outside the cached span.
    #[inline]
    pub fn states(&self, ids: &[i32], t: f64, out: &mut [[f64; 6]]) -> bool {
        self.states_rel(ids, t, 0.0, out)
    }

    /// [`PerturberCache::states`] at `jd_ref + dt` without rounding the sum to a Julian date.
    #[inline]
    pub fn states_rel(&self, ids: &[i32], jd_ref: f64, dt: f64, out: &mut [[f64; 6]]) -> bool {
        // Time since the start of the cache, formed without passing through an absolute JD.
        let rel = (jd_ref - self.t_start) + dt;
        if !(rel >= 0.0 && rel <= self.t_end - self.t_start) || ids != self.ids.as_slice() || out.len() < ids.len() {
            return false;
        }
        for (fit, o) in self.fits.iter().zip(out.iter_mut()) {
            let (k, x) = local_coordinate_rel(rel, fit.interval, fit.n_intervals);
            let block = &fit.coeffs[k * 6 * NCOEF..(k + 1) * 6 * NCOEF];
            // Chebyshev polynomials at x (a short scalar recurrence), then six independent dot
            // products. This has much less serial dependency than Clenshaw's recurrence.
            let mut tj = [0.0f64; NCOEF];
            tj[0] = 1.0;
            tj[1] = x;
            let x2 = 2.0 * x;
            for j in 2..NCOEF {
                tj[j] = x2 * tj[j - 1] - tj[j - 2];
            }
            let mut acc = [0.0f64; 6];
            for (j, cs) in block.chunks_exact(6).enumerate() {
                for c in 0..6 {
                    acc[c] += cs[c] * tj[j];
                }
            }
            *o = acc;
        }
        true
    }
}

/// Interval index and local coordinate in [-1, 1] of epoch `t`.
#[inline]
fn local_coordinate(t: f64, t_start: f64, interval: f64, n_intervals: usize) -> (usize, f64) {
    local_coordinate_rel(t - t_start, interval, n_intervals)
}

/// Interval index and local coordinate of a time `rel` days after the start of the fits.
#[inline]
fn local_coordinate_rel(rel: f64, interval: f64, n_intervals: usize) -> (usize, f64) {
    let u = rel / interval;
    let k = (u.floor() as isize).clamp(0, n_intervals as isize - 1) as usize;
    (k, 2.0 * (u - k as f64) - 1.0)
}

/// Evaluate piecewise Chebyshev fits laid out as `[interval][coefficient][component]`.
#[inline]
fn cheb_eval(coeffs: &[f64], t_start: f64, interval: f64, n_intervals: usize, t: f64, out: &mut [f64; 6]) {
    let (k, x) = local_coordinate(t, t_start, interval, n_intervals);
    let x2 = 2.0 * x;
    let block = &coeffs[k * 6 * NCOEF..(k + 1) * 6 * NCOEF];
    // Clenshaw recurrence, six components at once.
    let mut b1 = [0.0f64; 6];
    let mut b2 = [0.0f64; 6];
    for j in (1..NCOEF).rev() {
        let cs = &block[j * 6..j * 6 + 6];
        for c in 0..6 {
            let b0 = cs[c] + x2 * b1[c] - b2[c];
            b2[c] = b1[c];
            b1[c] = b0;
        }
    }
    for c in 0..6 {
        out[c] = block[c] + x * b1[c] - b2[c];
    }
}
