//! The orbit determination pipeline: layup's `do_fit` and `_orbitfit`.

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::orbfit::bk::bk_iod;
use crate::orbfit::herget::herget_iod;
use crate::orbfit::gauss::{gauss_astrometry, gauss_target_days, select_triplet};
use crate::orbfit::bk_fit::fit_orbit_bk;
use crate::orbfit::lm::fit_orbit;
use crate::orbfit::residuals::residuals;
use crate::orbfit::{Astrometry, Engine, FitFlag, FitOptions, NongravAuto, OrbitFit};
use crate::spice::SpiceKernel;

/// Smallest root Gauss's method keeps (AU).
const GAUSS_MIN_ROOT: f64 = 1e-4;
/// Physical bounds on the barycentric distance of an initial orbit (AU).
const MIN_R: f64 = 0.05;
const MAX_R: f64 = 1000.0;
/// Initial orbits passing closer than this to an observer (AU, in straight-line motion) are
/// fitted only if nothing else converges.
const CLOSE_APPROACH: f64 = 0.1;
/// Largest plausible hyperbolic excess speed (AU/day): 200 km/s.
const MAX_EXCESS_SPEED: f64 = 200.0 * 86400.0 / 149597870.7;

/// Split the detections into arcs separated by more than `gap` days and order the arcs for
/// fitting: the longest first, then repeatedly the one nearest in time to those already
/// included. Returns indices into `epochs`.
pub fn build_sequence(epochs: &[f64], gap: f64) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..epochs.len()).collect();
    order.sort_by(|&a, &b| epochs[a].total_cmp(&epochs[b]));
    let mut chunks: Vec<Vec<usize>> = Vec::new();
    for (k, &i) in order.iter().enumerate() {
        if k == 0 || epochs[i] - epochs[order[k - 1]] > gap {
            chunks.push(Vec::new());
        }
        chunks.last_mut().unwrap().push(i);
    }
    if chunks.is_empty() {
        return chunks;
    }
    let span = |c: &Vec<usize>| epochs[*c.last().unwrap()] - epochs[c[0]];
    let mut start = 0;
    for (k, c) in chunks.iter().enumerate() {
        if span(c) > span(&chunks[start]) {
            start = k;
        }
    }
    // Gap between two arcs.
    let distance = |a: &Vec<usize>, b: &Vec<usize>| {
        let (a0, a1) = (epochs[a[0]], epochs[*a.last().unwrap()]);
        let (b0, b1) = (epochs[b[0]], epochs[*b.last().unwrap()]);
        (a0 - b1).abs().min((b0 - a1).abs())
    };
    let mut seq = vec![chunks.remove(start)];
    while !chunks.is_empty() {
        let mut best = (f64::INFINITY, 0);
        for target in &seq {
            // nearest remaining arc to this one (first on ties)
            let mut nearest = (f64::INFINITY, usize::MAX);
            for (k, c) in chunks.iter().enumerate() {
                let d = distance(target, c);
                if d < nearest.0 {
                    nearest = (d, k);
                }
            }
            if nearest.0 < best.0 {
                best = nearest;
            }
        }
        seq.push(chunks.remove(best.1));
    }
    seq
}

/// Initial orbits from Gauss's method on a triplet drawn from `primary` (layup's `gauss_iod`).
fn gauss_candidates(astrometry: &Astrometry, primary: &[usize]) -> Vec<(f64, [f64; 6])> {
    let primary: Vec<usize> = primary.iter().copied().filter(|&i| astrometry.has_optical(i)).collect();
    if primary.len() < 3 {
        return Vec::new();
    }
    let (i, j, k) = select_triplet(&astrometry.epoch, &primary, gauss_target_days())
        .unwrap_or((primary[0], primary[primary.len() / 2], primary[primary.len() - 1]));
    gauss_astrometry(astrometry, i, j, k, GAUSS_MIN_ROOT)
}

/// Closest approach of a straight-line trajectory to the observers (AU).
fn min_observer_distance(astrometry: &Astrometry, epoch: f64, s: &[f64; 6]) -> f64 {
    let mut min_d2 = f64::INFINITY;
    for i in 0..astrometry.len() {
        let dt = astrometry.epoch[i] - epoch;
        let o = astrometry.observer[i];
        let d2 = (0..3).map(|k| (s[k] + s[k + 3] * dt - o[k]).powi(2)).sum::<f64>();
        if d2 < min_d2 {
            min_d2 = d2;
        }
    }
    min_d2.sqrt()
}

fn r2(s: &[f64; 6]) -> f64 {
    s[0] * s[0] + s[1] * s[1] + s[2] * s[2]
}

/// numpy's default (linear) percentile of unsorted data.
fn percentile(values: &mut [f64], q: f64) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    let pos = (values.len() - 1) as f64 * q / 100.0;
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(values.len() - 1);
    let t = pos - lo as f64;
    let (a, b) = (values[lo], values[hi]);
    if t >= 0.5 { b - (b - a) * (1.0 - t) } else { a + (b - a) * t }
}

/// Drop candidates that miss most detections by more than `opts.prefilter_sigma` (layup's
/// `filter_candidates_by_residual`). Never drops the best-fitting candidate.
fn prefilter(astrometry: &Astrometry, candidates: Vec<(f64, [f64; 6])>, kernel: &SpiceKernel, opts: &FitOptions) -> Vec<(f64, [f64; 6])> {
    let physical: Vec<(f64, [f64; 6])> =
        candidates.iter().copied().filter(|(_, s)| r2(s) > 0.0 && (MIN_R..=MAX_R).contains(&r2(s).sqrt())).collect();
    if physical.is_empty() {
        return candidates;
    }
    if astrometry.len() < 4 {
        return physical;
    }
    let mut keep = vec![false; physical.len()];
    let mut best: Option<(f64, usize)> = None;
    for (k, (epoch, s)) in physical.iter().enumerate() {
        if min_observer_distance(astrometry, *epoch, s) < CLOSE_APPROACH {
            keep[k] = true;
            continue;
        }
        let Ok(res) = residuals(astrometry, *epoch, s, &[0.0; 3], kernel, opts, false) else { continue };
        // Non-optical detections contribute zero, as in layup's residuals_at_state.
        let per = res.per_detection(astrometry.len());
        let mut sig: Vec<f64> = (0..astrometry.len())
            .map(|i| if astrometry.has_optical(i) { per[i][0].hypot(per[i][1]) / astrometry.sigma_ra[i].max(astrometry.sigma_dec[i]) } else { 0.0 })
            .collect();
        if sig.is_empty() {
            continue;
        }
        let metric = percentile(&mut sig, 80.0);
        if best.map_or(true, |(b, _)| metric < b) {
            best = Some((metric, k));
        }
        if metric <= opts.prefilter_sigma {
            keep[k] = true;
        }
    }
    let mut out: Vec<(f64, [f64; 6])> = physical.iter().zip(&keep).filter(|(_, &k)| k).map(|(c, _)| *c).collect();
    if let Some((_, k)) = best {
        if !keep[k] {
            out.push(physical[k]);
        }
    }
    if out.is_empty() { physical } else { out }
}

/// The converged fit with the smallest chi-square, preferring those farther than
/// `opts.min_distance` from the barycenter.
fn pick_best(fits: &[OrbitFit], opts: &FitOptions) -> Option<OrbitFit> {
    let converged: Vec<&OrbitFit> = fits.iter().filter(|f| f.converged()).collect();
    if converged.is_empty() {
        return None;
    }
    let sane: Vec<&OrbitFit> = converged.iter().copied().filter(|f| r2(&f.state) > opts.min_distance.powi(2)).collect();
    let pool = if sane.is_empty() { converged } else { sane };
    let mut best = pool[0];
    for f in &pool[1..] {
        if f.chi2 < best.chi2 {
            best = f;
        }
    }
    Some(best.clone())
}

/// layup's `_is_valid_data`: at least three detections, all values finite, none before 1801.
fn valid(a: &Astrometry) -> bool {
    const JD_1801: f64 = 2451545.0 - 6279962400.0 / 86400.0;
    a.len() >= 3
        && (0..a.len()).all(|i| {
            a.epoch[i] >= JD_1801 && (a.has_optical(i) || a.has_delay(i) || a.has_doppler(i)) && a.observer[i].iter().all(|v| v.is_finite())
        })
}

fn implausible(s: &[f64; 6]) -> bool {
    let r = r2(s).sqrt();
    if !s.iter().all(|v| v.is_finite()) || !(r > 0.0) {
        return false;
    }
    let energy = 0.5 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]) - GRAVITATIONAL_CONSTANT / r;
    energy > 0.0 && (2.0 * energy).sqrt() > MAX_EXCESS_SPEED
}

/// One gravity-only fit with `opts.engine` (layup's `_run_fit`). The BK engine ignores the
/// iteration budget and always allows `opts.max_iter`, as layup's does its fixed 100.
fn run_fit(a: &Astrometry, epoch: f64, state: &[f64; 6], kernel: &SpiceKernel, opts: &FitOptions, iters: usize) -> OrbitFit {
    match opts.engine {
        Engine::Cartesian => fit_orbit(a, epoch, state, &[0.0; 3], kernel, opts, iters),
        Engine::BkNative => fit_orbit_bk(a, epoch, state, kernel, opts, opts.max_iter),
    }
}

/// Gauss (or Herget) IOD, candidate screening, and differential correction (layup's `do_fit`), for a
/// gravity-only orbit.
fn gravity_fit(astrometry: &Astrometry, kernel: &SpiceKernel, opts: &FitOptions) -> OrbitFit {
    let seq = build_sequence(&astrometry.epoch, opts.arc_gap);
    if seq.is_empty() {
        return OrbitFit::failed(FitFlag::NoSolution);
    }
    let bk_fallback = opts.bk_fallback && !opts.herget;
    let mut solns = if opts.herget {
        let optical: Vec<usize> = seq[0].iter().copied().filter(|&i| astrometry.has_optical(i)).collect();
        herget_iod(astrometry, &optical, kernel, opts).into_iter().collect()
    } else {
        gauss_candidates(astrometry, &seq[0])
    };
    if solns.is_empty() && !bk_fallback {
        return OrbitFit::failed(FitFlag::NoSolution);
    }
    if solns.len() > 1 {
        solns = prefilter(astrometry, solns, kernel, opts);
    }

    // Fit every candidate on the primary arc with a small iteration budget; fall back to the
    // candidates passing close to an observer, then to the full budget.
    let primary = astrometry.subset(&seq[0]);
    let (safe, deferred): (Vec<_>, Vec<_>) =
        solns.iter().copied().partition(|(e, s)| min_observer_distance(astrometry, *e, s) >= CLOSE_APPROACH);
    let screen = |cands: &[(f64, [f64; 6])], iters: usize| -> Vec<OrbitFit> {
        cands.iter().map(|(e, s)| run_fit(&primary, *e, s, kernel, opts, iters)).collect()
    };
    let mut fits = screen(&safe, opts.screen_iter);
    let mut best = pick_best(&fits, opts);
    if best.is_none() && !deferred.is_empty() {
        fits.extend(screen(&deferred, opts.screen_iter));
        best = pick_best(&fits, opts);
    }
    if best.is_none() {
        fits = screen(&solns, opts.max_iter);
        best = pick_best(&fits, opts);
    }
    // Every Gauss candidate failed (or there were none): seed from the Bernstein-Khushalani
    // linear fit to the primary arc, at its middle detection.
    let primary_optical = primary.subset(&primary.optical());
    if best.is_none() && bk_fallback && primary_optical.len() >= 3 {
        let epoch = primary_optical.epoch[primary_optical.len() / 2];
        if let Some(seed) = bk_iod(&primary_optical, epoch) {
            let fit = run_fit(&primary, epoch, &seed, kernel, opts, opts.max_iter);
            if fit.converged() {
                best = Some(fit.clone());
            }
            fits.push(fit);
        }
    }
    let Some(primary_fit) = best else {
        if fits.is_empty() {
            return OrbitFit::failed(FitFlag::NoSolution);
        }
        // Report the least-bad attempt (first smallest chi-square, as Python's min()).
        let least_bad = fits.iter().fold(None::<&OrbitFit>, |b, f| match b {
            Some(b) if !(f.chi2 < b.chi2) => Some(b),
            _ => Some(f),
        });
        let mut out = least_bad.cloned().unwrap_or_else(|| OrbitFit::failed(FitFlag::NoSolution));
        out.flag = FitFlag::NoRootConverged;
        return out;
    };

    // All detections, starting from the primary-arc orbit.
    let fit = run_fit(astrometry, primary_fit.epoch, &primary_fit.state, kernel, opts, opts.max_iter);
    if fit.converged() {
        return fit;
    }
    // Otherwise add one arc at a time, starting again from the first candidate (or the
    // primary-arc orbit when it came from the Bernstein-Khushalani seed).
    let (mut epoch, mut state) = solns.first().copied().unwrap_or((primary_fit.epoch, primary_fit.state));
    let mut idx: Vec<usize> = Vec::new();
    let mut fit = fit;
    for arc in &seq {
        idx.extend(arc);
        fit = run_fit(&astrometry.subset(&idx), epoch, &state, kernel, opts, opts.max_iter);
        if !fit.converged() {
            fit.flag = FitFlag::BuildupFailed;
            return fit;
        }
        epoch = fit.epoch;
        state = fit.state;
    }
    // The last fit used every detection, in arc order: put the residuals back in input order.
    let mut resid = vec![[f64::NAN; 6]; fit.residuals.len()];
    for (k, &i) in idx.iter().enumerate() {
        resid[i] = fit.residuals[k];
    }
    fit.residuals = resid;
    fit
}

// ---- Robust mode (not in layup) -------------------------------------------------------------

/// A fit with a usable orbit, even if its chi-square is too large.
fn has_orbit(f: &OrbitFit) -> bool {
    matches!(f.flag, FitFlag::Converged | FitFlag::Chi2TooLarge) && f.state.iter().all(|v| v.is_finite())
}

/// Each detection's normalized residual: the rms over its measurements of residual / sigma.
fn normalized_residuals(a: &Astrometry, f: &OrbitFit, kernel: &SpiceKernel, opts: &FitOptions) -> Option<Vec<f64>> {
    let res = residuals(a, f.epoch, &f.state, &f.nongrav, kernel, opts, false).ok()?;
    let mut sum = vec![0.0; a.len()];
    let mut count = vec![0usize; a.len()];
    for ((&i, &k), &r) in res.row_detection.iter().zip(&res.row_kind).zip(&res.resid) {
        let z = r / a.sigma(i, k);
        sum[i] += z * z;
        count[i] += 1;
    }
    // NaN (a failed light-time solution, say) compares false below, so it counts as an outlier.
    Some((0..a.len()).map(|i| if count[i] > 0 { (sum[i] / count[i] as f64).sqrt() } else { f64::NAN }).collect())
}

/// Iterative outlier rejection starting from `start`, a fit to the detections `start_used`:
/// refit to the detections within `opts.outlier_sigma` of the current orbit until that set stops
/// changing. Rejected detections are re-evaluated every round. Returns the last fit that worked,
/// the detections it used, and whether the rejection settled: the used detections are exactly
/// those within `opts.outlier_sigma` of the final orbit (or, if the set still flips after the
/// last round, the fit's chi-square is within that bound on average).
fn reject_outliers(a: &Astrometry, start: OrbitFit, start_used: Vec<bool>, kernel: &SpiceKernel, opts: &FitOptions) -> (OrbitFit, Vec<bool>, bool) {
    const MAX_ROUNDS: usize = 10;
    let (mut fit, mut used) = (start, start_used);
    let mut refitted = false;
    for _ in 0..MAX_ROUNDS {
        let Some(z) = normalized_residuals(a, &fit, kernel, opts) else { return (fit, used, false) };
        let keep: Vec<bool> = z.iter().map(|&v| v <= opts.outlier_sigma).collect();
        if refitted && keep == used {
            return (fit, used, true);
        }
        let idx: Vec<usize> = (0..a.len()).filter(|&i| keep[i]).collect();
        if idx.len() < 3 {
            return (fit, used, false);
        }
        let refit = fit_orbit(&a.subset(&idx), fit.epoch, &fit.state, &fit.nongrav, kernel, opts, opts.max_iter);
        if !has_orbit(&refit) {
            return (fit, used, false);
        }
        fit = refit;
        used = keep;
        refitted = true;
    }
    let settled = refitted && fit.ndof > 0 && fit.chi2 <= opts.outlier_sigma.powi(2) * fit.ndof as f64;
    (fit, used, settled)
}

/// Orbits grown from short windows: the detections are split into windows of at most
/// `opts.seed_window` days, and the `opts.seed_tries` windows with the most detections are tried
/// in turn (most first). layup's pipeline (with `iod_opts`) gives each window an initial orbit,
/// which is then widened step by step to all detections, with outlier rejection at each step.
/// Returns the first orbit that reaches all detections.
fn grow_from_window(a: &Astrometry, kernel: &SpiceKernel, opts: &FitOptions, iod_opts: &FitOptions) -> Option<OrbitFit> {
    const MIN_WINDOW: usize = 6;
    let order = a.time_order();
    let t: Vec<f64> = order.iter().map(|&i| a.epoch[i]).collect();
    let n = t.len();
    let mut label = vec![0usize; n];
    let (mut w, mut start) = (0, t[0]);
    for k in 0..n {
        if t[k] - start > opts.seed_window {
            w += 1;
            start = t[k];
        }
        label[k] = w;
    }
    let mut count = vec![0usize; w + 1];
    for &l in &label {
        count[l] += 1;
    }
    let mut ranked: Vec<usize> = (0..=w).collect();
    ranked.sort_by(|&x, &y| count[y].cmp(&count[x])); // stable: earlier windows win ties

    'windows: for &win in ranked.iter().take(opts.seed_tries) {
        if count[win] < MIN_WINDOW {
            break;
        }
        let ks: Vec<usize> = (0..n).filter(|&k| label[k] == win).collect();
        let idx: Vec<usize> = ks.iter().map(|&k| order[k]).collect();
        let seed = gravity_fit(&a.subset(&idx), kernel, iod_opts);
        if !seed.converged() {
            continue;
        }
        let (mut fit, mut t0, mut t1) = (seed, t[ks[0]], t[*ks.last().unwrap()]);
        while t0 > t[0] || t1 < t[n - 1] {
            let span = (t1 - t0).max(10.0);
            t0 -= span;
            t1 += span;
            let idx: Vec<usize> = (0..a.len()).filter(|&i| a.epoch[i] >= t0 && a.epoch[i] <= t1).collect();
            let sub = a.subset(&idx);
            // New detections in, outliers and all; rejection then sorts them out.
            let refit = fit_orbit(&sub, fit.epoch, &fit.state, &[0.0; 3], kernel, opts, opts.max_iter);
            let start = if has_orbit(&refit) { refit } else { fit };
            let (f, _, settled) = reject_outliers(&sub, start, vec![true; sub.len()], kernel, opts);
            if !settled {
                continue 'windows;
            }
            fit = f;
        }
        return Some(fit);
    }
    None
}

/// Robust gravity-only fit. Stages, until one settles (see [`reject_outliers`]):
/// 1. layup's pipeline (or a fit from `initial`), then outlier rejection;
/// 2. orbits grown from short windows (see [`grow_from_window`]);
/// 3. both again with layup's chi-square test on initial orbits lifted, since outliers in the
///    arc used for the initial orbit fail every candidate otherwise.
/// If nothing settles, layup's own result is returned.
fn robust_gravity_fit(a: &Astrometry, initial: Option<(f64, [f64; 6])>, kernel: &SpiceKernel, opts: &FitOptions) -> (OrbitFit, Vec<bool>) {
    let all = vec![true; a.len()];
    let good = |f: &OrbitFit, used: &[bool], settled: bool| settled && has_orbit(f) && used.iter().filter(|&&u| u).count() >= 3;
    let relaxed = FitOptions { chi2_threshold: f64::INFINITY, ..opts.clone() };
    let mut plain = None;
    for iod_opts in [opts, &relaxed] {
        let first = match initial {
            Some((epoch, state)) => fit_orbit(a, epoch, &state, &[0.0; 3], kernel, iod_opts, opts.max_iter),
            None => gravity_fit(a, kernel, iod_opts),
        };
        // A fit that is merely pulled by outliers (chi-square too large, or a build-up that
        // failed on a bad arc) still starts the rejection.
        if has_orbit(&first) || (first.flag == FitFlag::BuildupFailed && first.state.iter().all(|v| v.is_finite())) {
            let (f, used, settled) = reject_outliers(a, first.clone(), all.clone(), kernel, opts);
            if good(&f, &used, settled) {
                return (f, used);
            }
        }
        if plain.is_none() {
            plain = Some(first); // layup's own result, to report if nothing settles
        }
        if initial.is_none() {
            if let Some(seed) = grow_from_window(a, kernel, opts, iod_opts) {
                let full = fit_orbit(a, seed.epoch, &seed.state, &[0.0; 3], kernel, opts, opts.max_iter);
                let start = if has_orbit(&full) { full } else { seed };
                let (f, used, settled) = reject_outliers(a, start, all.clone(), kernel, opts);
                if good(&f, &used, settled) {
                    return (f, used);
                }
            }
        }
    }
    (plain.unwrap(), all)
}

/// layup's model ladder for `fit_nongrav="auto"`: every single parameter before any pair, A2
/// first (issue #544 there). Within a tier, the largest chi-square drop wins.
const NONGRAV_LADDER: [&[[bool; 3]]; 3] = [
    &[[false, true, false], [true, false, false], [false, false, true]],
    &[[true, true, false]],
    &[[true, true, true]],
];

/// Automatic choice of the non-gravitational model (layup's `_select_nongrav_auto`), starting
/// from the converged gravity-only fit `gravity`. If its chi-square per degree of freedom is at
/// most `thresholds.accept_reduced_chi2`, it is returned as is. Otherwise the models A2, A1, A3,
/// then A1+A2, then A1+A2+A3 are fitted (tier by tier), and the first tier with a warranted
/// model gives the answer: the model with the smallest chi-square among those that converge
/// (flag 0), lower chi-square by more than `delta_chi2_per_param` per added parameter, and have
/// every added parameter above `nsigma` of its uncertainty. With none warranted, `gravity` is
/// returned.
pub fn select_nongrav(astrometry: &Astrometry, gravity: &OrbitFit, kernel: &SpiceKernel, opts: &FitOptions, thresholds: &NongravAuto) -> OrbitFit {
    if gravity.ndof <= 0 || gravity.chi2 / gravity.ndof as f64 <= thresholds.accept_reduced_chi2 {
        return gravity.clone();
    }
    for tier in NONGRAV_LADDER {
        let mut best: Option<OrbitFit> = None;
        for &mask in tier {
            let o = FitOptions { fit_nongrav: mask, nongrav_auto: None, nongrav_per_arc: false, ..opts.clone() };
            let f = fit_orbit(astrometry, gravity.epoch, &gravity.state, &[0.0; 3], kernel, &o, opts.max_iter);
            if f.flag != FitFlag::Converged {
                continue;
            }
            let k = mask.iter().filter(|&&b| b).count() as f64;
            if gravity.chi2 - f.chi2 <= thresholds.delta_chi2_per_param * k {
                continue;
            }
            let sigma = f.nongrav_sigma();
            if !(0..3).filter(|&j| mask[j]).all(|j| f.nongrav[j].abs() > thresholds.nsigma * sigma[j]) {
                continue;
            }
            if best.as_ref().map_or(true, |b| f.chi2 < b.chi2) {
                best = Some(f);
            }
        }
        if let Some(b) = best {
            return b;
        }
    }
    gravity.clone()
}

/// Determine an orbit from `astrometry`: layup's orbit fit.
///
/// Without an `initial` orbit (`(epoch, state, nongrav)`: TDB Julian date, barycentric J2000
/// state, and A1, A2, A3 to start a non-gravitational fit from), candidate orbits come from
/// Gauss's method on a triplet from the longest arc (see [`build_sequence`]); candidates that
/// badly miss the detections are discarded, the rest are fitted to the longest arc, and the best
/// converged one is fitted to all detections (or, failing that, to one more arc at a time).
/// With an `initial` orbit, only the final differential correction runs.
///
/// Non-gravitational parameters selected in `opts.fit_nongrav` are then fitted jointly,
/// starting from the gravity-only orbit; if that joint fit fails, the gravity-only orbit is
/// returned (with `fit_nongrav` all false).
///
/// With `opts.robust` (not in layup), detections more than `opts.outlier_sigma` from the orbit
/// are left out and re-evaluated each round until the set stops changing (for the gravity-only
/// and the non-gravitational fit alike). If layup's pipeline gives no orbit that way, the fit
/// starts over from the short window of detections that gives one and widens it to all
/// detections; failing that, both are retried with layup's chi-square test on initial orbits
/// lifted. `OrbitFit::used` marks the detections in the final fit, and `OrbitFit::residuals`
/// covers every detection, rejected ones included. If nothing works, the result is layup's.
pub fn determine_orbit(astrometry: &Astrometry, initial: Option<(f64, [f64; 6], [f64; 3])>, kernel: &SpiceKernel, opts: &FitOptions) -> OrbitFit {
    let mut fit = determine_orbit_inner(astrometry, initial, kernel, opts);
    // A failed fit may carry an attempt on a subset (the primary arc): keep per-detection
    // outputs aligned with the input.
    let n = astrometry.len();
    if fit.residuals.len() != n {
        fit.residuals = vec![[f64::NAN; 6]; n];
    }
    if fit.used.len() != n {
        fit.used = Vec::new();
    }
    fit
}

fn determine_orbit_inner(astrometry: &Astrometry, initial: Option<(f64, [f64; 6], [f64; 3])>, kernel: &SpiceKernel, opts: &FitOptions) -> OrbitFit {
    if !valid(astrometry) {
        return OrbitFit::failed(FitFlag::NotAttempted);
    }
    let mut gravity_opts = FitOptions { fit_nongrav: [false; 3], nongrav_auto: None, ..opts.clone() };
    // As layup: non-gravitational fitting runs entirely on the Cartesian engine.
    if opts.fit_nongrav.iter().any(|&b| b) || opts.nongrav_auto.is_some() {
        gravity_opts.engine = Engine::Cartesian;
    }
    if opts.robust {
        return robust_determine_orbit(astrometry, initial, kernel, opts, &gravity_opts);
    }
    let mut fit = match initial {
        Some((epoch, state, nongrav)) => {
            // As layup: the gravity-only fit ignores the initial non-gravitational parameters,
            // which then seed the joint fit.
            let mut f = run_fit(astrometry, epoch, &state, kernel, &gravity_opts, opts.max_iter);
            if opts.fit_nongrav.iter().any(|&b| b) {
                f.nongrav = nongrav;
            }
            f
        }
        None => gravity_fit(astrometry, kernel, &gravity_opts),
    };

    if let Some(thresholds) = opts.nongrav_auto {
        if fit.converged() {
            fit = select_nongrav(astrometry, &fit, kernel, opts, &thresholds);
        }
    } else if fit.converged() && opts.fit_nongrav.iter().any(|&b| b) {
        let ng = fit_orbit(astrometry, fit.epoch, &fit.state, &fit.nongrav, kernel, opts, opts.max_iter);
        if ng.converged() {
            fit = ng;
        }
    }
    if fit.converged() && implausible(&fit.state) {
        fit.flag = FitFlag::ImplausibleOrbit;
    }
    fit
}

fn robust_determine_orbit(a: &Astrometry, initial: Option<(f64, [f64; 6], [f64; 3])>, kernel: &SpiceKernel, opts: &FitOptions, gravity_opts: &FitOptions) -> OrbitFit {
    let (mut fit, mut used) = robust_gravity_fit(a, initial.map(|(e, s, _)| (e, s)), kernel, gravity_opts);
    // With automatic selection, the model is chosen on the detections that survived the
    // gravity-only fit, then fitted robustly like an explicit one.
    let mut opts = opts.clone();
    if let Some(thresholds) = opts.nongrav_auto.take() {
        opts.fit_nongrav = [false; 3];
        if fit.converged() {
            let idx: Vec<usize> = (0..a.len()).filter(|&i| used[i]).collect();
            let chosen = select_nongrav(&a.subset(&idx), &fit, kernel, &opts, &thresholds);
            if chosen.npar > 6 {
                opts.fit_nongrav = chosen.fit_nongrav;
            }
        }
    }
    let opts = &opts;
    if fit.converged() && opts.fit_nongrav.iter().any(|&b| b) {
        // The gravity-only rejection can drop detections that only the non-gravitational model
        // explains (precise radar of a Yarkovsky drifter, say), so start the joint fit both from
        // all detections and from the gravity-only set, and keep whichever settles on more.
        let start_ng = initial.map(|(_, _, ng)| ng).unwrap_or([0.0; 3]);
        let mut best: Option<(OrbitFit, Vec<bool>)> = None;
        let masks = [vec![true; a.len()], used.clone()];
        for mask in masks.iter() {
            let idx: Vec<usize> = (0..a.len()).filter(|&i| mask[i]).collect();
            let ng = fit_orbit(&a.subset(&idx), fit.epoch, &fit.state, &start_ng, kernel, opts, opts.max_iter);
            if !has_orbit(&ng) {
                continue;
            }
            let (f, u, settled) = reject_outliers(a, ng, mask.clone(), kernel, opts);
            if !(settled && f.converged()) {
                continue;
            }
            let n_used = u.iter().filter(|&&x| x).count();
            let better = match &best {
                None => true,
                Some((b, bu)) => {
                    let nb = bu.iter().filter(|&&x| x).count();
                    n_used > nb || (n_used == nb && f.chi2 < b.chi2)
                }
            };
            if better {
                best = Some((f, u));
            }
            if mask == &used {
                break; // both starts coincide when nothing was rejected
            }
        }
        if let Some((f, u)) = best {
            fit = f;
            used = u;
        }
    }
    if has_orbit(&fit) {
        // Residuals of every detection, rejected ones included, at the final orbit.
        if let Ok(res) = residuals(a, fit.epoch, &fit.state, &fit.nongrav, kernel, opts, false) {
            fit.residuals = res.per_detection(a.len());
        }
        fit.used = used;
    }
    if fit.converged() && implausible(&fit.state) {
        fit.flag = FitFlag::ImplausibleOrbit;
    }
    fit
}
