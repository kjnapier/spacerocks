//! Matching detections against a catalog of orbits.

use std::collections::HashMap;

use nalgebra::Vector3;
use rayon::prelude::*;

use crate::batch::BatchOptions;
use crate::checker::mpc::runoff_from_u;
use crate::checker::Catalog;
use crate::constants::SPEED_OF_LIGHT;
use crate::orbfit::predict::predict_with_states;
use crate::orbfit::{Astrometry, FitOptions};
use crate::spacerock::hg_magnitude;
use crate::transforms::universal_kepler_step;
use crate::{Origin, SpiceKernel};

type BoxError = Box<dyn std::error::Error + Send + Sync>;

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);
const DEGREE: f64 = std::f64::consts::PI / 180.0;
/// Days per decade (the MPC's runoff is per decade).
const DECADE: f64 = 3652.5;

/// Settings for [`check`].
#[derive(Debug, Clone)]
pub struct CheckOptions {
    /// A detection is consistent with an orbit when their Mahalanobis distance (on the sky, with
    /// the covariances of both) is at most this.
    pub nsigma: f64,
    /// Also report objects predicted within this angle (radians) of a detection even when they
    /// are not consistent with it, as MPChecker's search radius. 0 reports consistent pairs only.
    pub radius: f64,
    /// Two-body motion is used for the coarse search up to this many days from a reference
    /// state (the orbit's epoch, the catalog's snapshot, or a state integrated for the purpose).
    pub max_age: f64,
    /// Allowance (radians) for the error of two-body motion in the coarse search:
    /// `margin + growth * dt^2`, dt in days from the reference state, with `growth` =
    /// `margin_growth` for orbits with perihelion inside 1.3 AU and `margin_growth_distant`
    /// beyond. Over MPCORB, two-body motion from an N-body state drifts (as seen from the Earth)
    /// at most 5" in 5 days, 190" in 20 and 580" in 30 for perihelia inside 1.3 AU, and 0.2",
    /// 0.7" and 1.7" beyond; the defaults cover these with room to spare.
    pub margin: f64,
    pub margin_growth: f64,
    pub margin_growth_distant: f64,
    /// Objects whose predicted 1-sigma positional uncertainty exceeds this (radians) at a
    /// detection are only reported within `radius` of it (they are consistent with anything
    /// nearby).
    pub max_uncertainty: f64,
    /// Added in quadrature to every prediction's uncertainty (radians, per axis): what neither
    /// the orbit's covariance nor the detection's uncertainty covers, such as the rounding of
    /// catalog elements and force-model differences with the catalog's orbit computer.
    pub floor: f64,
    /// For orbits without a covariance (MPCORB), the uncertainty across the direction of motion
    /// as a fraction of the along-track uncertainty given by the U parameter.
    pub cross_track: f64,
    /// Detection epochs within this many days are searched together.
    pub window: f64,
    /// IAS15 tolerance of the N-body predictions.
    pub epsilon: f64,
    /// Catalog objects processed per block (bounds the memory of integrated states).
    pub block: usize,
    /// How catalog objects are integrated to the reference epochs the coarse search needs.
    pub batch: BatchOptions,
    pub parallel: bool,
}

impl Default for CheckOptions {
    fn default() -> Self {
        CheckOptions {
            nsigma: 3.0,
            radius: 60.0 * ARCSEC,
            max_age: 30.0,
            margin: 30.0 * ARCSEC,
            margin_growth: 1.0 * ARCSEC,
            margin_growth_distant: 0.02 * ARCSEC,
            max_uncertainty: 600.0 * ARCSEC,
            floor: 0.3 * ARCSEC,
            cross_track: 0.1,
            window: 1.0,
            epsilon: 1e-9,
            block: 100_000,
            batch: BatchOptions::default(),
            parallel: true,
        }
    }
}

/// A (detection, object) pair: the object's prediction at the detection, and how the two
/// compare.
#[derive(Debug, Clone, Copy)]
pub struct Match {
    /// Index of the detection.
    pub detection: usize,
    /// Index of the object in the catalog.
    pub object: usize,
    /// Predicted astrometric (light-time corrected) RA and Dec (radians).
    pub ra: f64,
    pub dec: f64,
    /// Observed minus predicted on the tangent plane at the prediction, along (RA cos Dec, Dec)
    /// (radians).
    pub offset: [f64; 2],
    /// Angle between the detection and the prediction (radians).
    pub separation: f64,
    /// The prediction's covariance (radians²) along (RA cos Dec, Dec).
    pub covariance: [[f64; 2]; 2],
    /// Mahalanobis distance of the offset, with the prediction's and the detection's
    /// covariances.
    pub distance: f64,
    /// `distance <= nsigma`.
    pub consistent: bool,
    /// Log of the Gaussian probability density of the offset (per arcsec²), with the
    /// prediction's and the detection's covariances: ranks a precise orbit that fits above a
    /// vague one that merely does not exclude the detection.
    pub log_likelihood: f64,
    /// Predicted rates dRA/dt and dDec/dt (radians/day; NaN without observer velocities).
    pub ra_rate: f64,
    pub dec_rate: f64,
    /// Observer–object distance and Sun–object distance at emission (AU).
    pub delta: f64,
    pub r_helio: f64,
    /// Predicted V magnitude (H, G system; NaN without H).
    pub mag: f64,
    /// Whether the uncertainty came from a covariance (true) or the MPC U parameter (false).
    pub from_covariance: bool,
}

impl Match {
    /// The prediction's 1-sigma error ellipse: semi-major and semi-minor axes (radians) and the
    /// position angle of the major axis (degrees, North through East).
    pub fn ellipse(&self) -> (f64, f64, f64) {
        crate::orbfit::Prediction { epoch: 0.0, ra: self.ra, dec: self.dec, delta: self.delta, cov: self.covariance }.ellipse()
    }
}

/// Result of [`check`].
#[derive(Debug, Clone, Default)]
pub struct Checked {
    /// Pairs that are consistent, or closer than [`CheckOptions::radius`], sorted by detection,
    /// then consistent pairs first, then by decreasing likelihood.
    pub matches: Vec<Match>,
    /// Objects whose prediction failed, with the reason.
    pub failed: Vec<(usize, String)>,
    /// Number of (detection, object) pairs that passed the coarse search.
    pub candidates: usize,
}

fn unit(ra: f64, dec: f64) -> Vector3<f64> {
    Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin())
}

fn radec(u: &Vector3<f64>) -> (f64, f64) {
    (u.y.atan2(u.x).rem_euclid(std::f64::consts::TAU), (u.z / u.norm()).clamp(-1.0, 1.0).asin())
}

fn angle(a: &Vector3<f64>, b: &Vector3<f64>) -> f64 {
    // Robust for small angles.
    2.0 * ((a - b).norm() / 2.0).clamp(0.0, 1.0).asin()
}

/// Two-body motion (heliocentric) by `dt` days, universal variables.
#[inline]
fn kepler(p: &Vector3<f64>, v: &Vector3<f64>, mu: f64, dt: f64) -> Option<(Vector3<f64>, Vector3<f64>)> {
    if dt == 0.0 {
        return Some((*p, *v));
    }
    let (p1, v1) = universal_kepler_step(p, v, mu, dt).ok()?;
    (p1.iter().all(|x| x.is_finite())).then_some((p1, v1))
}

/// Detections on a grid of cells on the sky, for radius queries.
struct SkyIndex {
    width: f64,
    bands: Vec<(usize, HashMap<usize, Vec<usize>>)>,
    all: Vec<usize>,
    units: Vec<Vector3<f64>>,
}

impl SkyIndex {
    fn new(items: &[usize], units: &[Vector3<f64>], width: f64) -> SkyIndex {
        let nb = (std::f64::consts::PI / width).ceil() as usize;
        let mut bands: Vec<(usize, HashMap<usize, Vec<usize>>)> = (0..nb)
            .map(|b| {
                let lo = -std::f64::consts::FRAC_PI_2 + b as f64 * width;
                let edge = lo.abs().min((lo + width).abs());
                let nra = ((std::f64::consts::TAU * edge.cos() / width).floor() as usize).max(1);
                (nra, HashMap::new())
            })
            .collect();
        for &i in items {
            let (ra, dec) = radec(&units[i]);
            let b = (((dec + std::f64::consts::FRAC_PI_2) / width) as usize).min(nb - 1);
            let nra = bands[b].0;
            let c = ((ra / std::f64::consts::TAU * nra as f64) as usize).min(nra - 1);
            bands[b].1.entry(c).or_default().push(i);
        }
        SkyIndex { width, bands, all: items.to_vec(), units: units.to_vec() }
    }

    /// Calls `f` for every item within `radius` of `u` (and possibly some a little farther).
    fn query(&self, u: &Vector3<f64>, radius: f64, mut f: impl FnMut(usize)) {
        if radius >= 20.0 * DEGREE || self.all.len() <= 16 {
            for &i in &self.all {
                if angle(u, &self.units[i]) <= radius {
                    f(i);
                }
            }
            return;
        }
        let (ra, dec) = radec(u);
        let nb = self.bands.len();
        let half = std::f64::consts::FRAC_PI_2;
        let b0 = (((dec - radius + half) / self.width).floor().max(0.0)) as usize;
        let b1 = ((((dec + radius + half) / self.width).floor()) as usize).min(nb - 1);
        for b in b0..=b1 {
            let (nra, cells) = &self.bands[b];
            let top = (dec.abs() + radius).min(half);
            let dra = if top >= half - 1e-9 { std::f64::consts::PI } else { (radius.sin() / top.cos()).min(1.0).asin() };
            let nra = *nra;
            if dra >= std::f64::consts::PI / 2.0 || nra <= 3 {
                for list in cells.values() {
                    for &i in list {
                        if angle(u, &self.units[i]) <= radius {
                            f(i);
                        }
                    }
                }
                continue;
            }
            let cw = std::f64::consts::TAU / nra as f64;
            let c0 = ((ra - dra) / cw).floor() as i64;
            let c1 = ((ra + dra) / cw).floor() as i64;
            for c in c0..=c1 {
                let c = c.rem_euclid(nra as i64) as usize;
                if let Some(list) = cells.get(&c) {
                    for &i in list {
                        if angle(u, &self.units[i]) <= radius {
                            f(i);
                        }
                    }
                }
            }
        }
    }
}

/// Detections sharing an epoch and an observer.
struct Group {
    t: f64,
    obs: Vector3<f64>,
    sun: [f64; 6],
    dets: Vec<usize>,
    /// Largest detection uncertainty in the group (radians).
    sigma: f64,
}

/// Groups within [`CheckOptions::window`] days, searched together.
struct Window {
    t: f64,
    groups: Vec<usize>,
    /// Observer position at the group nearest `t`, and the largest distance of another group's
    /// observer from it.
    obs: Vector3<f64>,
    obs_spread: f64,
    half_span: f64,
    index: SkyIndex,
    sigma: f64,
}

/// Windows whose epochs are within 2 * max_age days, sharing reference states.
struct Cluster {
    center: f64,
    lo: f64,
    hi: f64,
    windows: Vec<usize>,
    /// As for [`Window`], over the whole cluster.
    obs: Vector3<f64>,
    obs_spread: f64,
    index: SkyIndex,
    sigma: f64,
}

/// Perihelion distance of a heliocentric state (any conic).
fn perihelion(p: &Vector3<f64>, v: &Vector3<f64>, mu: f64) -> f64 {
    let r = p.norm();
    let h2 = p.cross(v).norm_squared();
    let e = ((v.norm_squared() - mu / r) * p - p.dot(v) * v).norm() / mu;
    h2 / (mu * (1.0 + e))
}

fn has_covariance(cat: &Catalog, i: usize) -> bool {
    let c = &cat.covariance[i];
    let n = (c.len() as f64).sqrt().round() as usize;
    n >= 6 && n * n == c.len()
}

/// A rough upper bound on the 1-sigma positional uncertainty (AU) of object `i` at `t`, for the
/// coarse search; `r_h` is its heliocentric distance there. With a covariance: the position's
/// plus four times the velocity's times the time from the epoch (the along-track drift a
/// velocity error causes). Without: the U parameter's along-track uncertainty.
fn sigma_hint(cat: &Catalog, i: usize, t: f64, r_h: f64) -> f64 {
    let c = &cat.covariance[i];
    if has_covariance(cat, i) {
        let dt = (t - cat.orbits.epochs[i]).abs();
        let n = (c.len() as f64).sqrt().round() as usize;
        let sp = (c[0] + c[n + 1] + c[2 * n + 2]).max(0.0).sqrt();
        let sv = (c[3 * n + 3] + c[4 * n + 4] + c[5 * n + 5]).max(0.0).sqrt();
        sp + 4.0 * sv * dt
    } else {
        r_h * longitude_sigma(cat, i, t)
    }
}

/// Sky covariance (radians², along RA cos Dec and Dec, at the direction `rho_hat`, distance
/// `delta`) of object `i` without a covariance: the U parameter's in-orbit longitude
/// uncertainty at `t`, as a shift `r_h * sigma` along the orbit (direction of the heliocentric
/// velocity `v_h`), plus `cross` times it across.
#[allow(clippy::too_many_arguments)]
fn u_sky_covariance(cat: &Catalog, i: usize, t: f64, r_h: f64, v_h: &Vector3<f64>, rho_hat: &Vector3<f64>, delta: f64, cross: f64) -> [[f64; 2]; 2] {
    let s = r_h * longitude_sigma(cat, i, t) / delta;
    let (a_vec, d_vec) = sky_basis(rho_hat);
    let vhat = v_h / v_h.norm();
    let w = [s * vhat.dot(&a_vec), s * vhat.dot(&d_vec)];
    let c2 = (cross * s) * (cross * s);
    [[w[0] * w[0] + c2, w[0] * w[1]], [w[1] * w[0], w[1] * w[1] + c2]]
}

/// Unit vectors along increasing RA and Dec at the direction `u`.
fn sky_basis(u: &Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let (ra, dec) = radec(u);
    (Vector3::new(-ra.sin(), ra.cos(), 0.0), Vector3::new(-dec.sin() * ra.cos(), -dec.sin() * ra.sin(), dec.cos()))
}

/// Whether the offset `x` lies within the ellipse `nsigma` sigma of `cov`, grown by `slack` on
/// each axis.
fn within(x: [f64; 2], cov: &[[f64; 2]; 2], nsigma: f64, slack: f64) -> bool {
    let [[a, b], [_, d]] = *cov;
    let tr = a + d;
    let disc = ((a - d).powi(2) + 4.0 * b * b).sqrt();
    let l1 = (0.5 * (tr + disc)).max(0.0);
    let l2 = (0.5 * (tr - disc)).max(0.0);
    // Eigenvector of l1.
    let (ex, ey) = if b.abs() > 1e-300 { (l1 - d, b) } else if a >= d { (1.0, 0.0) } else { (0.0, 1.0) };
    let norm = (ex * ex + ey * ey).sqrt();
    let (ex, ey) = (ex / norm, ey / norm);
    let p1 = x[0] * ex + x[1] * ey;
    let p2 = -x[0] * ey + x[1] * ex;
    let s1 = nsigma * l1.sqrt() + slack;
    let s2 = nsigma * l2.sqrt() + slack;
    (p1 / s1).powi(2) + (p2 / s2).powi(2) <= 1.0
}

/// In-orbit longitude uncertainty (radians) of object `i` at `t`, from the MPC uncertainty
/// parameter U (unknown counts as U = 9): the runoff per decade, growing linearly with the time
/// from the epoch of the elements (or from the last observation, when that is longer), plus a
/// tenth of it (a year's growth). Growth from the epoch also covers what the elements
/// themselves lose (their rounding to 7 digits, forces the catalog's orbit included and this
/// integration does not): real detections of (101955) Bennu and (99942) Apophis during Earth
/// approaches 15-20 years before the epoch are off by up to ~30".
fn longitude_sigma(cat: &Catalog, i: usize, t: f64) -> f64 {
    let u = cat.u[i];
    let runoff = runoff_from_u(if u.is_finite() { u } else { 9.0 });
    let mut dt = (t - cat.orbits.epochs[i]).abs();
    if cat.last_obs[i].is_finite() {
        dt = dt.max(t - cat.last_obs[i]);
    }
    runoff * (0.1 + dt / DECADE)
}

/// Check `detections` (epochs, RA/Dec, uncertainties and observers of [`Astrometry`]; RA
/// uncertainty along RA cos Dec) against every orbit in `catalog`. `correlation` is the RA–Dec
/// correlation of each detection's uncertainty (empty for none).
///
/// Two stages:
///
/// 1. **Coarse**, over the whole catalog: objects move by two-body motion (about the Sun) from
///    a reference state within [`CheckOptions::max_age`] days: the orbit's own epoch, the
///    catalog's snapshot, or a state integrated (N-body, all objects of the block at once) to
///    the middle of the detections' epochs. An object is a candidate for a detection when its
///    predicted position is within `max(radius, nsigma * sigma) + margin` of it, sigma
///    combining the detection's uncertainty and a generous bound on the object's.
/// 2. **Refined**, for each candidate object: an N-body prediction (ASSIST's force model,
///    light-time corrected, as [`crate::orbfit::predict`]) at each of its candidate
///    detections, with the covariance mapped through the variational equations; or, for orbits
///    without a covariance, the along-track uncertainty of the MPC's U parameter projected along
///    the object's motion. The detection is consistent with the orbit when the Mahalanobis
///    distance of their offset, `sqrt(r^T (C_pred + C_det)^-1 r)`, is at most `nsigma`.
pub fn check(catalog: &Catalog, detections: &Astrometry, correlation: &[f64], kernel: &SpiceKernel, opts: &CheckOptions) -> Result<Checked, BoxError> {
    let nd = detections.epoch.len();
    if detections.ra.len() != nd || detections.dec.len() != nd || detections.observer.len() != nd {
        return Err("detections: epoch, ra, dec and observer must have the same length".into());
    }
    if !correlation.is_empty() && correlation.len() != nd {
        return Err(format!("{} correlations for {} detections", correlation.len(), nd).into());
    }
    let sig = |v: &Vec<f64>, j: usize| v.get(j).copied().filter(|s| s.is_finite() && *s > 0.0).unwrap_or(crate::orbfit::DEFAULT_SIGMA);
    let det_units: Vec<Vector3<f64>> = (0..nd).map(|j| unit(detections.ra[j], detections.dec[j])).collect();
    let det_sigma: Vec<f64> = (0..nd).map(|j| sig(&detections.sigma_ra, j).max(sig(&detections.sigma_dec, j))).collect();

    // Groups (epoch, observer), in time order.
    let mut gmap: HashMap<(u64, [u64; 3]), usize> = HashMap::new();
    let mut groups: Vec<Group> = Vec::new();
    let mut order: Vec<usize> = (0..nd).filter(|&j| detections.epoch[j].is_finite() && detections.ra[j].is_finite() && detections.dec[j].is_finite()).collect();
    order.sort_by(|&a, &b| detections.epoch[a].total_cmp(&detections.epoch[b]));
    for &j in &order {
        let o = detections.observer[j];
        let key = (detections.epoch[j].to_bits(), [o[0].to_bits(), o[1].to_bits(), o[2].to_bits()]);
        let g = *gmap.entry(key).or_insert_with(|| {
            groups.push(Group { t: detections.epoch[j], obs: Vector3::from(o), sun: [0.0; 6], dets: Vec::new(), sigma: 0.0 });
            groups.len() - 1
        });
        groups[g].dets.push(j);
        groups[g].sigma = groups[g].sigma.max(det_sigma[j]);
    }
    for g in groups.iter_mut() {
        g.sun = kernel.state_au(10, 0, g.t)?;
    }
    if groups.is_empty() {
        return Ok(Checked::default());
    }

    // Windows, then clusters of windows.
    let mut windows: Vec<Window> = Vec::new();
    let mut start = 0;
    while start < groups.len() {
        let mut end = start + 1;
        while end < groups.len() && groups[end].t - groups[start].t <= opts.window {
            end += 1;
        }
        let ids: Vec<usize> = (start..end).collect();
        let t = 0.5 * (groups[start].t + groups[end - 1].t);
        let near = *ids.iter().min_by(|&&a, &&b| (groups[a].t - t).abs().total_cmp(&(groups[b].t - t).abs())).unwrap();
        let obs = groups[near].obs;
        let obs_spread = ids.iter().map(|&g| (groups[g].obs - obs).norm()).fold(0.0, f64::max);
        let dets: Vec<usize> = ids.iter().flat_map(|&g| groups[g].dets.iter().copied()).collect();
        let sigma = ids.iter().map(|&g| groups[g].sigma).fold(0.0, f64::max);
        windows.push(Window { t, half_span: 0.5 * (groups[end - 1].t - groups[start].t), groups: ids, obs, obs_spread, index: SkyIndex::new(&dets, &det_units, DEGREE), sigma });
        start = end;
    }
    let mut clusters: Vec<Cluster> = Vec::new();
    let mut start = 0;
    while start < windows.len() {
        let lo = groups[windows[start].groups[0]].t;
        let mut end = start + 1;
        while end < windows.len() && groups[*windows[end].groups.last().unwrap()].t - lo <= 2.0 * opts.max_age {
            end += 1;
        }
        let hi = groups[*windows[end - 1].groups.last().unwrap()].t;
        let center = 0.5 * (lo + hi);
        let ids: Vec<usize> = (start..end).flat_map(|w| windows[w].groups.iter().copied()).collect();
        let near = *ids.iter().min_by(|&&a, &&b| (groups[a].t - center).abs().total_cmp(&(groups[b].t - center).abs())).unwrap();
        let obs = groups[near].obs;
        let obs_spread = ids.iter().map(|&g| (groups[g].obs - obs).norm()).fold(0.0, f64::max);
        let dets: Vec<usize> = ids.iter().flat_map(|&g| groups[g].dets.iter().copied()).collect();
        let sigma = (start..end).map(|w| windows[w].sigma).fold(0.0, f64::max);
        clusters.push(Cluster { center, lo, hi, windows: (start..end).collect(), obs, obs_spread, index: SkyIndex::new(&dets, &det_units, DEGREE), sigma });
        start = end;
    }

    let mu = Origin::SUN.mu();
    let snap = catalog.has_snapshot();
    let n = catalog.len();
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    let all: Vec<usize> = (0..n).collect();
    for block in all.chunks(opts.block.max(1)) {
        // Reference state per (object, cluster): (epoch, state), or an index into `integrated`.
        let age = |t: f64, c: &Cluster| (c.lo - t).abs().max((c.hi - t).abs());
        let mut need: Vec<Vec<usize>> = vec![Vec::new(); clusters.len()];
        for &i in block {
            for (c, cl) in clusters.iter().enumerate() {
                let own = age(catalog.orbits.epochs[i], cl);
                let sn = if snap && catalog.snapshot[i][0].is_finite() { age(catalog.snapshot_epoch, cl) } else { f64::INFINITY };
                if own.min(sn) > opts.max_age {
                    need[c].push(i);
                }
            }
        }
        let need_c: Vec<usize> = (0..clusters.len()).filter(|&c| !need[c].is_empty()).collect();
        let mut need_objs: Vec<usize> = need_c.iter().flat_map(|&c| need[c].iter().copied()).collect();
        need_objs.sort_unstable();
        need_objs.dedup();
        let targets: Vec<f64> = need_c.iter().map(|&c| clusters[c].center).collect();
        let integrated = if need_objs.is_empty() { Vec::new() } else { catalog.states_at_nearest(&need_objs, &targets, kernel, &opts.batch) };
        let obj_row: HashMap<usize, usize> = need_objs.iter().enumerate().map(|(k, &i)| (i, k)).collect();
        let tgt_col: HashMap<usize, usize> = need_c.iter().enumerate().map(|(k, &c)| (c, k)).collect();
        let mt = targets.len();

        // The Sun at every reference epoch.
        let mut suns: HashMap<u64, [f64; 6]> = HashMap::new();
        let mut ref_epochs: Vec<f64> = block.iter().map(|&i| catalog.orbits.epochs[i]).collect();
        ref_epochs.extend(&targets);
        if snap {
            ref_epochs.push(catalog.snapshot_epoch);
        }
        for t in ref_epochs {
            if t.is_finite() && !suns.contains_key(&t.to_bits()) {
                if let Ok(s) = kernel.state_au(10, 0, t) {
                    suns.insert(t.to_bits(), s);
                }
            }
        }

        // Whether object `i` (heliocentric reference state `p_r`, `v_r` at `t_r`) could appear
        // near any detection of `index`, whose epochs are within `half` days of `t` and whose
        // observers are within `spread` AU of `obs`. The Sun's state at `t` is `sun`.
        #[allow(clippy::too_many_arguments)]
        let near = |i: usize, p_r: &Vector3<f64>, v_r: &Vector3<f64>, t_r: f64, growth: f64, t: f64, half: f64, obs: &Vector3<f64>, spread: f64, sun: &[f64; 6], index: &SkyIndex, sigma: f64| -> bool {
            let Some((p_w, v_w)) = kepler(p_r, v_r, mu, t - t_r) else { return false };
            let rho_w = p_w + Vector3::new(sun[0], sun[1], sun[2]) - obs;
            let delta_w = rho_w.norm();
            let t_far = if (t + half - catalog.orbits.epochs[i]).abs() > (t - half - catalog.orbits.epochs[i]).abs() { t + half } else { t - half };
            let sigma_obj = sigma_hint(catalog, i, t_far, p_w.norm()) / delta_w;
            let dt_max = (t - t_r).abs() + half;
            // Objects too uncertain to be matched are still reported within `radius`.
            let consistency = if sigma_obj > opts.max_uncertainty { 0.0 } else { opts.nsigma * (sigma_obj * sigma_obj + sigma * sigma + opts.floor * opts.floor).sqrt() };
            let gate = opts.radius.max(consistency)
                + opts.margin
                + growth * dt_max * dt_max;
            // How far the object can appear to move (its motion, the observers' spread, the
            // Sun's motion, and light time), as an angle.
            let r = p_w.norm();
            let shift = v_w.norm() * half + 0.5 * mu / (r * r) * half * half + spread + 0.02 * half + delta_w / SPEED_OF_LIGHT * v_w.norm();
            let bound = if shift < 0.9 * delta_w { (shift / (delta_w - shift)).atan() * 1.05 } else { std::f64::consts::PI };
            let mut hit = false;
            index.query(&(rho_w / delta_w), bound + gate, |_| hit = true);
            hit
        };

        let coarse = |&i: &usize| -> Vec<(usize, usize)> {
            let mut out = Vec::new();
            for (c, cl) in clusters.iter().enumerate() {
                // The reference state.
                let (t_r, s_r) = if let (Some(&k), Some(&col)) = (obj_row.get(&i), tgt_col.get(&c)) {
                    (targets[col], integrated[k * mt + col])
                } else {
                    let own = age(catalog.orbits.epochs[i], cl);
                    let sn = if snap && catalog.snapshot[i][0].is_finite() { age(catalog.snapshot_epoch, cl) } else { f64::INFINITY };
                    if sn < own { (catalog.snapshot_epoch, catalog.snapshot[i]) } else { (catalog.orbits.epochs[i], catalog.orbits.states[i]) }
                };
                if !s_r.iter().all(|x| x.is_finite()) {
                    continue;
                }
                let Some(sun_r) = suns.get(&t_r.to_bits()) else { continue };
                let p_r = Vector3::new(s_r[0] - sun_r[0], s_r[1] - sun_r[1], s_r[2] - sun_r[2]);
                let v_r = Vector3::new(s_r[3] - sun_r[3], s_r[4] - sun_r[4], s_r[5] - sun_r[5]);
                let growth = if perihelion(&p_r, &v_r, mu) < 1.3 { opts.margin_growth } else { opts.margin_growth_distant };
                if cl.windows.len() > 1 {
                    let sun_c = &groups[windows[cl.windows[0]].groups[0]].sun;
                    let half = (cl.center - cl.lo).max(cl.hi - cl.center);
                    if !near(i, &p_r, &v_r, t_r, growth, cl.center, half, &cl.obs, cl.obs_spread, sun_c, &cl.index, cl.sigma) {
                        continue;
                    }
                }
                for &w in &cl.windows {
                    let win = &windows[w];
                    let sun_w = &groups[win.groups[0]].sun;
                    if !near(i, &p_r, &v_r, t_r, growth, win.t, win.half_span, &win.obs, win.obs_spread, sun_w, &win.index, win.sigma) {
                        continue;
                    }
                    // Each group of the window, exactly.
                    for &g in &win.groups {
                        let grp = &groups[g];
                        let Some((p, v)) = kepler(&p_r, &v_r, mu, grp.t - t_r) else { continue };
                        let sun = Vector3::new(grp.sun[0], grp.sun[1], grp.sun[2]);
                        let vb = v + Vector3::new(grp.sun[3], grp.sun[4], grp.sun[5]);
                        let rho = p + sun - grp.obs;
                        let lt = rho.norm() / SPEED_OF_LIGHT;
                        let rho = rho - vb * lt;
                        let delta = rho.norm();
                        let u = rho / delta;
                        let dt = (grp.t - t_r).abs();
                        let slack = opts.margin + growth * dt * dt;
                        let cov_obj = has_covariance(catalog, i);
                        let s_obj = sigma_hint(catalog, i, grp.t, p.norm()) / delta;
                        let too_uncertain = s_obj > opts.max_uncertainty;
                        let (a_vec, d_vec) = sky_basis(&u);
                        let c_u = if cov_obj { [[0.0; 2]; 2] } else { u_sky_covariance(catalog, i, grp.t, p.norm(), &v, &u, delta, opts.cross_track) };
                        for &j in &grp.dets {
                            let sep = angle(&u, &det_units[j]);
                            if sep <= opts.radius + slack {
                                out.push((j, i));
                                continue;
                            }
                            if too_uncertain {
                                continue;
                            }
                            let s2 = det_sigma[j] * det_sigma[j] + opts.floor * opts.floor;
                            let ok = if cov_obj {
                                sep <= opts.nsigma * (s_obj * s_obj + s2).sqrt() + slack
                            } else {
                                let dj = det_units[j];
                                let x = [dj.dot(&a_vec), dj.dot(&d_vec)];
                                within(x, &[[c_u[0][0] + s2, c_u[0][1]], [c_u[1][0], c_u[1][1] + s2]], opts.nsigma, slack)
                            };
                            if ok {
                                out.push((j, i));
                            }
                        }
                    }
                }
            }
            out
        };
        let found: Vec<Vec<(usize, usize)>> = if opts.parallel { block.par_iter().map(coarse).collect() } else { block.iter().map(coarse).collect() };
        pairs.extend(found.into_iter().flatten());
    }
    let candidates = pairs.len();

    // Refinement, one integration per object.
    let mut by_obj: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(j, i) in &pairs {
        by_obj.entry(i).or_default().push(j);
    }
    let mut objs: Vec<(usize, Vec<usize>)> = by_obj.into_iter().collect();
    objs.sort_by_key(|(i, _)| *i);
    let fit_opts = FitOptions { epsilon: opts.epsilon, parallel: false, ..FitOptions::default() };
    let rho_of = |j: usize| correlation.get(j).copied().filter(|r| r.is_finite()).unwrap_or(0.0);

    let refine = |(i, dets): &(usize, Vec<usize>)| -> Result<Vec<Match>, (usize, String)> {
        let i = *i;
        let has_cov = has_covariance(catalog, i);
        let mut dets = dets.clone();
        dets.sort_by(|&a, &b| detections.epoch[a].total_cmp(&detections.epoch[b]));
        dets.dedup();
        // Start from the orbit's own state (where its covariance is defined), or, without a
        // covariance, from the snapshot if that is closer to the detections.
        let t_mid = 0.5 * (detections.epoch[dets[0]] + detections.epoch[*dets.last().unwrap()]);
        let (t0, s0) = if !has_cov && snap && catalog.snapshot[i][0].is_finite() && (catalog.snapshot_epoch - t_mid).abs() < (catalog.orbits.epochs[i] - t_mid).abs() {
            (catalog.snapshot_epoch, catalog.snapshot[i])
        } else {
            (catalog.orbits.epochs[i], catalog.orbits.states[i])
        };
        let fit = catalog.orbit_fit(i, s0, t0);
        let epochs: Vec<f64> = dets.iter().map(|&j| detections.epoch[j]).collect();
        let observers: Vec<[f64; 3]> = dets.iter().map(|&j| detections.observer[j]).collect();
        let preds = predict_with_states(&fit, &epochs, &observers, kernel, &fit_opts, has_cov).map_err(|e| (i, e))?;
        let mut out = Vec::new();
        for (k, &j) in dets.iter().enumerate() {
            let (p, obj) = &preds[k];
            let rho_hat = unit(p.ra, p.dec);
            let a_vec = Vector3::new(-p.ra.sin(), p.ra.cos(), 0.0);
            let d_vec = Vector3::new(-p.dec.sin() * p.ra.cos(), -p.dec.sin() * p.ra.sin(), p.dec.cos());
            let t_emit = p.epoch - p.delta / SPEED_OF_LIGHT;
            let sun = kernel.state_au(10, 0, t_emit).map_err(|e| (i, e.to_string()))?;
            let pos = Vector3::new(obj[0], obj[1], obj[2]);
            let vel = Vector3::new(obj[3], obj[4], obj[5]);
            let helio = pos - Vector3::new(sun[0], sun[1], sun[2]);
            let r_helio = helio.norm();
            let mut cov = p.cov;
            if !has_cov {
                // Along-track: a shift of the object along its (heliocentric) orbit.
                let vh = vel - Vector3::new(sun[3], sun[4], sun[5]);
                let c = u_sky_covariance(catalog, i, p.epoch, r_helio, &vh, &rho_hat, p.delta, opts.cross_track);
                for u in 0..2 {
                    for v in 0..2 {
                        cov[u][v] += c[u][v];
                    }
                }
            }
            cov[0][0] += opts.floor * opts.floor;
            cov[1][1] += opts.floor * opts.floor;
            let d = det_units[j];
            let cosang = d.dot(&rho_hat);
            let offset = [d.dot(&a_vec) / cosang, d.dot(&d_vec) / cosang];
            let separation = angle(&d, &rho_hat);
            let (sa, sd) = (sig(&detections.sigma_ra, j), sig(&detections.sigma_dec, j));
            let rho = rho_of(j);
            let s = [[cov[0][0] + sa * sa, cov[0][1] + rho * sa * sd], [cov[1][0] + rho * sa * sd, cov[1][1] + sd * sd]];
            let detm = s[0][0] * s[1][1] - s[0][1] * s[1][0];
            let d2 = (s[1][1] * offset[0] * offset[0] - (s[0][1] + s[1][0]) * offset[0] * offset[1] + s[0][0] * offset[1] * offset[1]) / detm;
            let distance = if cosang > 0.0 && detm > 0.0 { d2.max(0.0).sqrt() } else { f64::INFINITY };
            let log_likelihood = -0.5 * distance * distance - (std::f64::consts::TAU * detm.max(0.0).sqrt() / (ARCSEC * ARCSEC)).ln();
            let consistent = distance <= opts.nsigma;
            let sigma_major = crate::orbfit::Prediction { epoch: p.epoch, ra: p.ra, dec: p.dec, delta: p.delta, cov }.ellipse().0;
            if !(separation <= opts.radius || (consistent && sigma_major <= opts.max_uncertainty)) {
                continue;
            }
            let (ra_rate, dec_rate) = match detections.observer_velocity.get(j) {
                Some(vo) => {
                    let rel = vel - Vector3::from(*vo);
                    let cd = p.dec.cos();
                    (rel.dot(&a_vec) / (p.delta * cd), rel.dot(&d_vec) / p.delta)
                }
                None => (f64::NAN, f64::NAN),
            };
            let obs = Vector3::from(detections.observer[j]);
            let phase = helio.angle(&(pos - obs));
            let mag = if catalog.h[i].is_finite() { hg_magnitude(catalog.h[i], catalog.g[i], r_helio, p.delta, phase) } else { f64::NAN };
            out.push(Match {
                detection: j,
                object: i,
                ra: p.ra,
                dec: p.dec,
                offset,
                separation,
                covariance: cov,
                distance,
                consistent,
                log_likelihood,
                ra_rate,
                dec_rate,
                delta: p.delta,
                r_helio,
                mag,
                from_covariance: has_cov,
            });
        }
        Ok(out)
    };
    let results: Vec<Result<Vec<Match>, (usize, String)>> = if opts.parallel { objs.par_iter().map(refine).collect() } else { objs.iter().map(refine).collect() };
    let mut out = Checked { candidates, ..Default::default() };
    for r in results {
        match r {
            Ok(m) => out.matches.extend(m),
            Err(e) => out.failed.push(e),
        }
    }
    out.matches.sort_by(|a, b| a.detection.cmp(&b.detection).then(b.consistent.cmp(&a.consistent)).then(b.log_likelihood.total_cmp(&a.log_likelihood)));
    Ok(out)
}
