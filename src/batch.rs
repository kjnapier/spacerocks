//! Batch propagation and ephemerides for many rocks and many epochs.
//!
//! The functions here are the fast path for population-scale work. Compared with calling
//! [`SpaceRock::propagate`] and [`SpaceRock::observe`] in a loop they
//!
//! * integrate many rocks in one IAS15 simulation, so the perturber ephemeris is evaluated once
//!   per step for the whole group instead of once per rock;
//! * integrate each group once through all requested epochs (sorted), reading states off the
//!   integrator's dense output instead of restarting the integration for every epoch;
//! * convert frames and origins with one ephemeris lookup per distinct epoch; and
//! * run groups in parallel.
//!
//! Rocks sharing a simulation share its adaptive step size, so rocks are grouped by epoch and
//! then by perihelion distance before being split into chunks of [`BatchOptions::chunk_size`]:
//! near-Earth objects end up with other near-Earth objects and do not force small steps on
//! TNOs. The shared step is never longer than any member's own would be, so in practice each
//! rock is integrated at least as accurately as alone; results differ from single-rock
//! integration at the level of the integrator's tolerance.

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;

use crate::assist::{PerturberCache, SpiceSimulation};
use crate::observing::{apparent, Apparent};
use crate::{Observer, Origin, ReferencePlane, SpaceRock, SpiceKernel, Time};

type BoxError = Box<dyn std::error::Error + Send + Sync>;

/// How rocks are moved between epochs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// IAS15 with the Sun, planets, Moon, Pluto and 16 massive asteroids from the kernel
    /// (the same force model as [`SpaceRock::propagate`]).
    NBody,
    /// Two-body (Keplerian) motion about each rock's own origin, as in
    /// [`SpaceRock::analytic_propagate`].
    TwoBody,
}

impl Method {
    pub fn from_str(s: &str) -> Result<Method, String> {
        match s.to_lowercase().as_str() {
            "nbody" | "n-body" | "ias15" => Ok(Method::NBody),
            "twobody" | "two-body" | "kepler" | "keplerian" | "analytic" => Ok(Method::TwoBody),
            _ => Err(format!("unknown method '{}' (expected 'nbody' or 'twobody')", s)),
        }
    }
}

/// Options for [`propagate_batch`] and [`ephemeris`].
#[derive(Debug, Clone)]
pub struct BatchOptions {
    /// Maximum number of rocks integrated together in one simulation. When running in
    /// parallel, smaller groups are made if needed to keep every thread busy (down to 8 rocks).
    pub chunk_size: usize,
    pub method: Method,
    /// Fit the perturber ephemeris once over the time span and share it between all chunks
    /// (see [`PerturberCache`]). Off by default: once rocks share simulations, the ephemeris is
    /// a small part of the cost and building the cache does not pay for itself.
    pub perturber_cache: bool,
    /// Run chunks in parallel (rayon).
    pub parallel: bool,
    /// Also return each rock's barycentric J2000 state at every epoch ([`ephemeris`] only).
    pub with_states: bool,
}

impl Default for BatchOptions {
    fn default() -> Self {
        BatchOptions { chunk_size: 64, method: Method::NBody, perturber_cache: false, parallel: true, with_states: false }
    }
}

/// Result of [`ephemeris`]: `n_rocks x n_epochs` arrays, stored row-major (rock-major).
#[derive(Debug, Clone)]
pub struct Ephemeris {
    pub n_rocks: usize,
    pub n_epochs: usize,
    /// TDB Julian dates of the epochs, in the order given.
    pub epochs: Vec<f64>,
    pub apparent: Vec<Apparent>,
    /// Barycentric J2000 states (AU, AU/day) at each epoch, if requested.
    pub states: Option<Vec<[f64; 6]>>,
}

impl Ephemeris {
    pub fn get(&self, rock: usize, epoch: usize) -> &Apparent {
        &self.apparent[rock * self.n_epochs + epoch]
    }
}

/// Matrix taking vectors in `plane` to J2000: the inverse of the plane's matrix, exactly as
/// `SpaceRock::change_reference_plane` computes it (not every plane matrix is orthonormal to
/// machine precision, so the transpose would not round-trip).
fn to_j2000(plane: &ReferencePlane) -> Matrix3<f64> {
    plane.get_rotation_matrix().try_inverse().unwrap_or_else(Matrix3::identity)
}

fn state6(p: &Vector3<f64>, v: &Vector3<f64>) -> [f64; 6] {
    [p.x, p.y, p.z, v.x, v.y, v.z]
}

fn pv(s: &[f64; 6]) -> (Vector3<f64>, Vector3<f64>) {
    (Vector3::new(s[0], s[1], s[2]), Vector3::new(s[3], s[4], s[5]))
}

/// Sun's barycentric J2000 state at each distinct TDB Julian date.
struct SunCache {
    map: HashMap<u64, [f64; 6]>,
}

impl SunCache {
    fn new(jds: impl Iterator<Item = f64>, kernel: &SpiceKernel) -> Result<SunCache, BoxError> {
        let mut map = HashMap::new();
        for jd in jds {
            if let std::collections::hash_map::Entry::Vacant(e) = map.entry(jd.to_bits()) {
                e.insert(kernel.state_au(10, 0, jd)?);
            }
        }
        Ok(SunCache { map })
    }

    fn get(&self, jd: f64) -> [f64; 6] {
        self.map[&jd.to_bits()]
    }
}

/// Barycentric J2000 state of a rock at its own epoch.
fn rock_to_ssb_j2000(rock: &SpaceRock, suns: &SunCache, kernel: &SpiceKernel) -> Result<[f64; 6], BoxError> {
    match &rock.origin {
        Origin::SSB | Origin::SUN => {
            let m = to_j2000(&rock.reference_plane);
            let mut p = m * rock.position;
            let mut v = m * rock.velocity;
            if rock.origin == Origin::SUN {
                let (sp, sv) = pv(&suns.get(rock.epoch.tdb().jd()));
                p += sp;
                v += sv;
            }
            Ok(state6(&p, &v))
        }
        Origin::Custom { .. } => {
            let mut r = rock.clone();
            r.change_reference_plane("J2000").map_err(|e| e.to_string())?;
            r.to_ssb(kernel).map_err(|e| format!("{}: {}", rock.name, e))?;
            Ok(state6(&r.position, &r.velocity))
        }
    }
}

/// Barycentric J2000 observer (position, velocity, Sun position) at its epoch.
fn observer_to_ssb_j2000(o: &Observer, kernel: &SpiceKernel) -> Result<(Vector3<f64>, Vector3<f64>, Vector3<f64>), BoxError> {
    let m = to_j2000(&o.reference_plane);
    let mut p = m * o.position;
    let mut v = m * o.velocity.unwrap_or_else(Vector3::zeros);
    let jd = o.epoch.tdb().jd();
    let (sun_p, sun_v) = pv(&kernel.state_au(10, 0, jd)?);
    match &o.origin {
        Origin::SSB => {}
        Origin::SUN => {
            p += sun_p;
            v += sun_v;
        }
        Origin::Custom { name, .. } => {
            return Err(format!("observers with a custom origin ('{}') are not supported; use SSB or SUN", name).into())
        }
    }
    Ok((p, v, sun_p))
}

/// Perihelion distance about the SSB (grouping key; any conic).
fn perihelion_key(s: &[f64; 6]) -> f64 {
    let mu = Origin::SSB.mu();
    let (p, v) = pv(s);
    let h2 = p.cross(&v).norm_squared();
    let r = p.norm();
    let e = ((v.norm_squared() - mu / r) * p - p.dot(&v) * v).norm() / mu;
    let q = h2 / (mu * (1.0 + e));
    if q.is_finite() { q } else { f64::MAX }
}

/// A group of rocks with a common starting epoch, integrated together.
struct Chunk {
    t0: f64,
    rocks: Vec<usize>,
}

fn make_chunks(t0s: &[f64], states: &[[f64; 6]], active: &[usize], chunk_size: usize, parallel: bool) -> Vec<Chunk> {
    let mut groups: HashMap<u64, Vec<usize>> = HashMap::new();
    for &i in active {
        groups.entry(t0s[i].to_bits()).or_default().push(i);
    }
    let mut keys: Vec<u64> = groups.keys().copied().collect();
    keys.sort_by(|a, b| f64::from_bits(*a).total_cmp(&f64::from_bits(*b)));
    let mut chunks = Vec::new();
    let chunk_size = chunk_size.max(1);
    let threads = if parallel { rayon::current_num_threads() } else { 1 };
    for k in keys {
        let mut idx = groups.remove(&k).unwrap();
        let chunk_size = chunk_size.min(idx.len().div_ceil(threads).max(8));
        idx.sort_by(|&a, &b| perihelion_key(&states[a]).total_cmp(&perihelion_key(&states[b])));
        for c in idx.chunks(chunk_size) {
            chunks.push(Chunk { t0: f64::from_bits(k), rocks: c.to_vec() });
        }
    }
    chunks
}

/// Integrate one chunk through the (unique, sorted) target epochs. Returns the states as
/// `[rock_in_chunk][target]`.
fn run_nbody_chunk(
    chunk: &Chunk,
    rocks: &[SpaceRock],
    states: &[[f64; 6]],
    targets: &[f64],
    kernel: &SpiceKernel,
    cache: Option<&Arc<PerturberCache>>,
) -> Result<Vec<[f64; 6]>, BoxError> {
    let n = chunk.rocks.len();
    let m = targets.len();
    let mut out = vec![[f64::NAN; 6]; n * m];

    // Forward and backward in time from t0 are separate integrations.
    let split = targets.partition_point(|&t| t < chunk.t0);
    let backward: Vec<usize> = (0..split).rev().collect();
    let forward: Vec<usize> = (split..m).collect();

    for order in [forward, backward] {
        if order.is_empty() {
            continue;
        }
        let epoch = Time::new(chunk.t0, "tdb", "jd").map_err(|e| e.to_string())?;
        let mut sim = SpiceSimulation::horizons(&epoch, kernel).map_err(|e| e.to_string())?;
        if let Some(c) = cache {
            sim.set_perturber_cache(Arc::clone(c));
        }
        for &i in &chunk.rocks {
            let (p, v) = pv(&states[i]);
            let rock = SpaceRock {
                name: String::new(),
                epoch: epoch.clone(),
                reference_plane: ReferencePlane::J2000,
                origin: Origin::SSB,
                position: p,
                velocity: v,
                // Physical properties carry non-gravitational parameters (and mass).
                properties: rocks[i].properties.clone(),
            };
            sim.add(rock).map_err(|e| e.to_string())?;
        }
        for &j in &order {
            sim.integrate_jd(targets[j], kernel).map_err(|e| e.to_string())?;
            for (k, p) in sim.state.particles.iter().enumerate() {
                out[k * m + j] = state6(&p.position, &p.velocity);
            }
        }
    }
    Ok(out)
}

fn unique_sorted(v: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut u: Vec<f64> = v.to_vec();
    u.sort_by(|a, b| a.total_cmp(b));
    u.dedup();
    let index = v.iter().map(|x| u.binary_search_by(|y| y.total_cmp(x)).unwrap()).collect();
    (u, index)
}

fn build_cache(
    t0s: &[f64],
    targets: &[f64],
    kernel: &SpiceKernel,
    opts: &BatchOptions,
    n_chunks: usize,
) -> Option<Arc<PerturberCache>> {
    if !opts.perturber_cache || opts.method != Method::NBody || n_chunks < 2 {
        return None;
    }
    let lo = t0s.iter().chain(targets).copied().fold(f64::INFINITY, f64::min);
    let hi = t0s.iter().chain(targets).copied().fold(f64::NEG_INFINITY, f64::max);
    if !(lo.is_finite() && hi.is_finite()) || hi <= lo {
        return None;
    }
    let ids = SpiceSimulation::horizons_body_ids().ok()?;
    // Integration steps may overshoot the last target, so cover some margin if the kernels allow.
    let margin = 0.05 * (hi - lo) + 60.0;
    PerturberCache::build(kernel, &ids, lo - margin, hi + margin)
        .or_else(|_| PerturberCache::build(kernel, &ids, lo, hi))
        .ok()
        .map(Arc::new)
}

/// Barycentric J2000 states of every rock at every target epoch (TDB Julian dates), row-major
/// `[rock][target]`, by `opts.method`. Each rock is integrated once through all targets.
pub fn states_at(
    rocks: &[SpaceRock],
    targets: &[f64],
    kernel: &SpiceKernel,
    opts: &BatchOptions,
) -> Result<Vec<[f64; 6]>, BoxError> {
    let n = rocks.len();
    let m = targets.len();
    let mut out = vec![[f64::NAN; 6]; n * m];
    if n == 0 || m == 0 {
        return Ok(out);
    }
    let t0s: Vec<f64> = rocks.iter().map(|r| r.epoch.tdb().jd()).collect();

    match opts.method {
        Method::TwoBody => {
            let suns = SunCache::new(targets.iter().copied(), kernel)?;
            let target_times: Vec<Time> = targets
                .iter()
                .map(|&t| Time::new(t, "tdb", "jd"))
                .collect::<Result<_, _>>()
                .map_err(|e| e.to_string())?;
            let work = |(i, row): (usize, &mut [[f64; 6]])| -> Result<(), BoxError> {
                let rock = &rocks[i];
                if let Origin::Custom { name, .. } = &rock.origin {
                    return Err(format!("{}: two-body batch propagation needs a SUN or SSB origin, not '{}'", rock.name, name).into());
                }
                let m_rot = to_j2000(&rock.reference_plane);
                for (j, t) in target_times.iter().enumerate() {
                    let r = rock.analytic_at(t).map_err(|e| format!("{}: {}", rock.name, e))?;
                    let mut p = m_rot * r.position;
                    let mut v = m_rot * r.velocity;
                    if rock.origin == Origin::SUN {
                        let (sp, sv) = pv(&suns.get(targets[j]));
                        p += sp;
                        v += sv;
                    }
                    row[j] = state6(&p, &v);
                }
                Ok(())
            };
            if opts.parallel {
                out.par_chunks_mut(m).enumerate().try_for_each(work)?;
            } else {
                out.chunks_mut(m).enumerate().try_for_each(work)?;
            }
        }
        Method::NBody => {
            let suns = SunCache::new(t0s.iter().copied(), kernel)?;
            let states: Vec<[f64; 6]> = if opts.parallel {
                rocks.par_iter().map(|r| rock_to_ssb_j2000(r, &suns, kernel)).collect::<Result<_, _>>()?
            } else {
                rocks.iter().map(|r| rock_to_ssb_j2000(r, &suns, kernel)).collect::<Result<_, _>>()?
            };
            let all: Vec<usize> = (0..n).collect();
            let chunks = make_chunks(&t0s, &states, &all, opts.chunk_size, opts.parallel);
            let cache = build_cache(&t0s, targets, kernel, opts, chunks.len());
            let run = |c: &Chunk| {
                run_nbody_chunk(c, rocks, &states, targets, kernel, cache.as_ref()).map_err(|e| -> BoxError {
                    let names: Vec<&str> = c.rocks.iter().take(3).map(|&i| rocks[i].name.as_str()).collect();
                    let more = if c.rocks.len() > 3 { format!(" and {} more", c.rocks.len() - 3) } else { String::new() };
                    format!("integrating {}{}: {}", names.join(", "), more, e).into()
                })
            };
            let results: Vec<Vec<[f64; 6]>> = if opts.parallel {
                chunks.par_iter().map(run).collect::<Result<_, _>>()?
            } else {
                chunks.iter().map(run).collect::<Result<_, _>>()?
            };
            for (c, r) in chunks.iter().zip(results) {
                for (k, &i) in c.rocks.iter().enumerate() {
                    out[i * m..(i + 1) * m].copy_from_slice(&r[k * m..(k + 1) * m]);
                }
            }
        }
    }
    Ok(out)
}

/// Propagate many rocks to one epoch.
///
/// Equivalent to calling [`SpaceRock::propagate`] (for [`Method::NBody`]) or
/// [`SpaceRock::analytic_propagate`] (for [`Method::TwoBody`]) on each rock, but much faster for
/// large collections. Each rock keeps its reference plane; SUN and SSB origins are kept, and
/// custom origins are returned as SSB (as `propagate` does). Rocks already at `epoch` are left
/// untouched.
pub fn propagate_batch(rocks: &mut [SpaceRock], epoch: &Time, kernel: &SpiceKernel, opts: &BatchOptions) -> Result<(), BoxError> {
    let t = epoch.tdb().jd();
    let active: Vec<usize> = (0..rocks.len()).filter(|&i| rocks[i].epoch.tdb().jd() != t).collect();
    if active.is_empty() {
        return Ok(());
    }

    if opts.method == Method::TwoBody {
        let work = |r: &mut SpaceRock| {
            if r.epoch.tdb().jd() == t {
                return Ok(());
            }
            r.analytic_propagate(epoch).map_err(|e| BoxError::from(format!("{}: {}", r.name, e)))
        };
        return if opts.parallel {
            rocks.par_iter_mut().try_for_each(work)
        } else {
            rocks.iter_mut().try_for_each(work)
        };
    }

    let subset: Vec<SpaceRock> = active.iter().map(|&i| rocks[i].clone()).collect();
    let states = states_at(&subset, &[t], kernel, opts)?;
    let (sun_p, sun_v) = pv(&kernel.state_au(10, 0, t)?);

    for (k, &i) in active.iter().enumerate() {
        let rock = &mut rocks[i];
        let (mut p, mut v) = pv(&states[k]);
        let origin = match rock.origin {
            Origin::SUN => {
                p -= sun_p;
                v -= sun_v;
                Origin::SUN
            }
            _ => Origin::SSB,
        };
        let m = rock.reference_plane.get_rotation_matrix();
        rock.position = m * p;
        rock.velocity = m * v;
        rock.origin = origin;
        rock.epoch = epoch.clone();
    }
    Ok(())
}

/// Ephemerides of many rocks at many epochs.
///
/// `observers` gives one observer per epoch (for example from
/// [`crate::Observatory::at`]); their epochs are the output epochs, in any order, and may be
/// before or after the rocks' epochs. Each rock is integrated once through all epochs.
/// RA/Dec are in the J2000 (ICRF) equator regardless of the rocks' or observers' reference
/// planes. The rocks themselves are not modified.
pub fn ephemeris(rocks: &[SpaceRock], observers: &[Observer], kernel: &SpiceKernel, opts: &BatchOptions) -> Result<Ephemeris, BoxError> {
    let n = rocks.len();
    let m = observers.len();
    let epochs: Vec<f64> = observers.iter().map(|o| o.epoch.tdb().jd()).collect();
    let (targets, target_index) = unique_sorted(&epochs);

    let obs: Vec<(Vector3<f64>, Vector3<f64>, Vector3<f64>)> =
        observers.iter().map(|o| observer_to_ssb_j2000(o, kernel)).collect::<Result<_, _>>()?;

    let states = states_at(rocks, &targets, kernel, opts)?;
    let mt = targets.len();

    let mut apps = vec![Apparent::NAN; n * m];
    let mut out_states = if opts.with_states { Some(vec![[f64::NAN; 6]; n * m]) } else { None };

    let fill = |(i, row): (usize, &mut [Apparent])| {
        let rock = &rocks[i];
        let (h, g) = match &rock.properties {
            Some(p) => (p.absolute_magnitude, p.gslope.unwrap_or(0.15)),
            None => (None, 0.15),
        };
        for j in 0..m {
            let (p, v) = pv(&states[i * mt + target_index[j]]);
            let (op, ov, sun) = &obs[j];
            row[j] = apparent(&p, &v, op, ov, sun, h, g);
        }
    };
    if opts.parallel {
        apps.par_chunks_mut(m.max(1)).enumerate().for_each(fill);
    } else {
        apps.chunks_mut(m.max(1)).enumerate().for_each(fill);
    }
    if let Some(s) = out_states.as_mut() {
        for i in 0..n {
            for j in 0..m {
                s[i * m + j] = states[i * mt + target_index[j]];
            }
        }
    }

    Ok(Ephemeris { n_rocks: n, n_epochs: m, epochs, apparent: apps, states: out_states })
}
