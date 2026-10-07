//! A catalog of orbits to check detections against.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::batch::{states_at, BatchOptions};
use crate::checker::mpc::{self, MpcElements, MpcOrb};
use crate::orbfit::{FitFlag, OrbitFit};
use crate::{Origin, Population, Properties, ReferencePlane, SpaceRock, SpiceKernel, Time};

type BoxError = Box<dyn std::error::Error + Send + Sync>;

/// Orbits, stored as columns (one entry per object).
///
/// Each orbit is a barycentric J2000 state (AU, AU/day) at a TDB Julian date, with, when known,
/// its covariance (row-major `npar x npar`: the state, then the fitted non-gravitational
/// parameters, as in [`OrbitFit`]). Orbits without a covariance (MPCORB) carry the MPC
/// uncertainty parameter `u` instead, which the check turns into an along-track uncertainty.
///
/// A catalog may also hold a *snapshot*: every orbit integrated (N-body) to one common epoch.
/// The check moves objects by two-body motion from whichever reference (the orbit's own epoch
/// or the snapshot) is closest in time, so a snapshot near the detections saves integrating the
/// whole catalog on every call. Snapshots are kept by [`Catalog::save`].
#[derive(Debug, Clone)]
pub struct Catalog {
    /// Names, epochs (TDB Julian dates) and barycentric J2000 states of the orbits; the other
    /// columns line up with these.
    pub orbits: Population,
    /// Non-gravitational parameters A1–A3 (AU/day²).
    pub nongrav: Vec<[f64; 3]>,
    /// Which of `nongrav` have rows in `covariance`.
    pub fit_nongrav: Vec<[bool; 3]>,
    /// Row-major covariance of the state and the fitted non-gravitational parameters; empty if
    /// unknown.
    pub covariance: Vec<Vec<f64>>,
    /// Absolute magnitude and slope parameter (H NaN if unknown).
    pub h: Vec<f64>,
    pub g: Vec<f64>,
    /// MPC uncertainty parameter U (0–9), NaN if unknown.
    pub u: Vec<f64>,
    /// Julian date of the last observation in the orbit (NaN if unknown): the uncertainty from
    /// U grows from there.
    pub last_obs: Vec<f64>,
    /// TDB Julian date of `snapshot` (NaN if there is none).
    pub snapshot_epoch: f64,
    /// Barycentric J2000 states at `snapshot_epoch` (empty if there is none; NaN rows where the
    /// integration failed).
    pub snapshot: Vec<[f64; 6]>,
}

impl Default for Catalog {
    fn default() -> Self {
        Catalog {
            orbits: Population::new(ReferencePlane::J2000, Origin::SSB),
            nongrav: Vec::new(),
            fit_nongrav: Vec::new(),
            covariance: Vec::new(),
            h: Vec::new(),
            g: Vec::new(),
            u: Vec::new(),
            last_obs: Vec::new(),
            snapshot_epoch: f64::NAN,
            snapshot: Vec::new(),
        }
    }
}

struct SunStates<'a> {
    kernel: &'a SpiceKernel,
    map: HashMap<u64, [f64; 6]>,
}

impl<'a> SunStates<'a> {
    fn new(kernel: &'a SpiceKernel) -> Self {
        SunStates { kernel, map: HashMap::new() }
    }

    fn get(&mut self, jd: f64) -> Result<[f64; 6], BoxError> {
        if let Some(s) = self.map.get(&jd.to_bits()) {
            return Ok(*s);
        }
        let s = self.kernel.state_au(10, 0, jd)?;
        self.map.insert(jd.to_bits(), s);
        Ok(s)
    }
}

fn add6(a: &[f64; 6], b: &[f64; 6]) -> [f64; 6] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4], a[5] + b[5]]
}

impl Catalog {
    pub fn len(&self) -> usize {
        self.orbits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.orbits.is_empty()
    }

    pub fn has_snapshot(&self) -> bool {
        self.snapshot.len() == self.len() && self.snapshot_epoch.is_finite()
    }

    #[allow(clippy::too_many_arguments)]
    fn push(&mut self, name: String, epoch: f64, state: [f64; 6], nongrav: [f64; 3], fit_nongrav: [bool; 3], covariance: Vec<f64>, h: f64, g: f64, u: f64, last_obs: f64) {
        self.orbits.names.push(name);
        self.orbits.epochs.push(epoch);
        self.orbits.states.push(state);
        self.orbits.properties.push(None);
        self.nongrav.push(nongrav);
        self.fit_nongrav.push(fit_nongrav);
        self.covariance.push(covariance);
        self.h.push(h);
        self.g.push(g);
        self.u.push(u);
        self.last_obs.push(last_obs);
    }

    /// Orbits from MPC elements (heliocentric ecliptic, TT), made barycentric J2000. Orbits
    /// whose elements do not give a state (e.g. e >= 1 with a > 0) are left out.
    pub fn from_mpc_elements(elements: &[MpcElements], kernel: &SpiceKernel) -> Result<Catalog, BoxError> {
        let mut cat = Catalog::default();
        let mut sun = SunStates::new(kernel);
        for el in elements {
            let Some(helio) = el.helio_j2000() else { continue };
            let t = el.epoch_tdb();
            let s = add6(&helio, &sun.get(t)?);
            cat.push(el.name.clone(), t, s, [0.0; 3], [false; 3], Vec::new(), el.h, el.g, el.u, el.last_obs);
        }
        Ok(cat)
    }

    /// The MPC's catalog of minor-planet orbits (MPCORB). `path` may be `mpcorb_extended.json`,
    /// `MPCORB.DAT`, or either gzipped; `None` uses `mpcorb_extended.json.gz` in
    /// [`mpc::mpc_dir`], downloading it first if it is missing and `download` is set. To get a
    /// newer catalog, call [`mpc::mpcorb_path`] with `update`.
    pub fn mpcorb(path: Option<&Path>, download: bool, kernel: &SpiceKernel) -> Result<Catalog, BoxError> {
        let path = match path {
            Some(p) => p.to_path_buf(),
            None => mpc::mpcorb_path(download, false)?,
        };
        // Parsing the JSON takes ~30 s; the states are cached next to it in binary
        // (`<file>.srcat`), used while it is newer than the file.
        let mut cache = path.clone().into_os_string();
        cache.push(".srcat");
        let cache = std::path::PathBuf::from(cache);
        let newer = |a: &Path, b: &Path| -> bool {
            match (std::fs::metadata(a).and_then(|m| m.modified()), std::fs::metadata(b).and_then(|m| m.modified())) {
                (Ok(ta), Ok(tb)) => ta >= tb,
                _ => false,
            }
        };
        if newer(&cache, &path) {
            if let Ok(c) = Catalog::load(&cache) {
                return Ok(c);
            }
        }
        let mut cat = Catalog::from_mpc_elements(&mpc::read_mpcorb(&path)?, kernel)?;
        // About 1% of the orbits have older epochs than the rest: a snapshot at the common epoch
        // saves integrating them on every check.
        if let Some(t) = cat.common_epoch() {
            cat.make_snapshot(t, kernel, &BatchOptions::default(), 100_000);
        }
        // Best effort: an unwritable directory just means no cache.
        let _ = cat.save(&cache);
        Ok(cat)
    }

    /// Orbits from MPC `mpc_orb` JSON records (with covariances).
    pub fn from_mpc_orbs(orbs: &[MpcOrb], kernel: &SpiceKernel) -> Result<Catalog, BoxError> {
        let mut cat = Catalog::default();
        let mut sun = SunStates::new(kernel);
        for o in orbs {
            let s = add6(&o.helio_j2000, &sun.get(o.epoch_tdb)?);
            // The Sun's state is exact, so the covariance of the barycentric state is the
            // heliocentric one.
            cat.push(o.name.clone(), o.epoch_tdb, s, o.nongrav, [false; 3], o.covariance.to_vec(), o.h, o.g, o.u, f64::NAN);
        }
        Ok(cat)
    }

    /// Orbits from fits (with their covariances, including fitted non-gravitational
    /// parameters). Fits without an orbit are left out; per-arc fits are not supported.
    pub fn from_fits(names: &[String], fits: &[OrbitFit], h: Option<&[f64]>) -> Result<Catalog, BoxError> {
        if names.len() != fits.len() {
            return Err(format!("{} names for {} fits", names.len(), fits.len()).into());
        }
        let mut cat = Catalog::default();
        for (k, (name, f)) in names.iter().zip(fits).enumerate() {
            if !f.state.iter().all(|x| x.is_finite()) || !f.epoch.is_finite() {
                continue;
            }
            if f.per_arc {
                return Err(format!("{}: per-arc fits are not supported in a catalog", name).into());
            }
            let cov = if f.covariance.len() == f.npar * f.npar && f.covariance.iter().all(|x| x.is_finite()) { f.covariance.clone() } else { Vec::new() };
            let hk = h.and_then(|h| h.get(k).copied()).unwrap_or(f64::NAN);
            cat.push(name.clone(), f.epoch, f.state, f.nongrav, f.fit_nongrav, cov, hk, 0.15, f64::NAN, f64::NAN);
        }
        Ok(cat)
    }

    /// Append another catalog's orbits (dropping both snapshots unless they are at the same
    /// epoch).
    pub fn extend(&mut self, other: &Catalog) {
        let keep = self.has_snapshot() && other.has_snapshot() && self.snapshot_epoch == other.snapshot_epoch;
        self.orbits.names.extend_from_slice(&other.orbits.names);
        self.orbits.epochs.extend_from_slice(&other.orbits.epochs);
        self.orbits.states.extend_from_slice(&other.orbits.states);
        self.orbits.properties.extend_from_slice(&other.orbits.properties);
        self.nongrav.extend_from_slice(&other.nongrav);
        self.fit_nongrav.extend_from_slice(&other.fit_nongrav);
        self.covariance.extend_from_slice(&other.covariance);
        self.h.extend_from_slice(&other.h);
        self.g.extend_from_slice(&other.g);
        self.u.extend_from_slice(&other.u);
        self.last_obs.extend_from_slice(&other.last_obs);
        if keep {
            self.snapshot.extend_from_slice(&other.snapshot);
        } else {
            self.snapshot.clear();
            self.snapshot_epoch = f64::NAN;
        }
    }

    /// The orbits at `indices`, in that order.
    pub fn select(&self, indices: &[usize]) -> Catalog {
        let pick = |v: &Vec<f64>| indices.iter().map(|&i| v[i]).collect::<Vec<_>>();
        Catalog {
            orbits: self.orbits.select(indices),
            nongrav: indices.iter().map(|&i| self.nongrav[i]).collect(),
            fit_nongrav: indices.iter().map(|&i| self.fit_nongrav[i]).collect(),
            covariance: indices.iter().map(|&i| self.covariance[i].clone()).collect(),
            h: pick(&self.h),
            g: pick(&self.g),
            u: pick(&self.u),
            last_obs: pick(&self.last_obs),
            snapshot_epoch: if self.has_snapshot() { self.snapshot_epoch } else { f64::NAN },
            snapshot: if self.has_snapshot() { indices.iter().map(|&i| self.snapshot[i]).collect() } else { Vec::new() },
        }
    }

    /// The most frequent epoch, if there are orbits.
    pub fn common_epoch(&self) -> Option<f64> {
        let mut counts: HashMap<u64, usize> = HashMap::new();
        for t in &self.orbits.epochs {
            *counts.entry(t.to_bits()).or_default() += 1;
        }
        counts.into_iter().max_by_key(|&(_, c)| c).map(|(t, _)| f64::from_bits(t))
    }

    /// Index of each name (the first, if a name repeats).
    pub fn index(&self) -> HashMap<&str, usize> {
        let mut m = HashMap::with_capacity(self.len());
        for (i, n) in self.orbits.names.iter().enumerate().rev() {
            m.insert(n.as_str(), i);
        }
        m
    }

    /// Object `i` as a barycentric J2000 SpaceRock at its epoch (with H, G and non-grav).
    pub fn rock(&self, i: usize) -> Result<SpaceRock, BoxError> {
        let t = Time::new(self.orbits.epochs[i], "tdb", "jd").map_err(|e| e.to_string())?;
        let s = self.orbits.states[i];
        let mut r = SpaceRock::from_xyz(&self.orbits.names[i], s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").map_err(|e| e.to_string())?;
        if self.nongrav[i].iter().any(|&a| a != 0.0) {
            r.set_nongrav(self.nongrav[i][0], self.nongrav[i][1], self.nongrav[i][2]);
        }
        if self.h[i].is_finite() {
            r.set_absolute_magnitude(self.h[i]);
            r.set_gslope(self.g[i]);
        }
        Ok(r)
    }

    /// Object `i` as an [`OrbitFit`] (converged, with its covariance or NaN).
    pub(crate) fn orbit_fit(&self, i: usize, state: [f64; 6], epoch: f64) -> OrbitFit {
        let cov = &self.covariance[i];
        let npar = (cov.len() as f64).sqrt().round() as usize;
        let (covariance, npar) = if npar >= 6 && npar * npar == cov.len() { (cov.clone(), npar) } else { (vec![f64::NAN; 36], 6) };
        let mut f = OrbitFit::failed(FitFlag::Converged);
        f.epoch = epoch;
        f.state = state;
        f.nongrav = self.nongrav[i];
        f.fit_nongrav = if npar > 6 { self.fit_nongrav[i] } else { [false; 3] };
        f.covariance = covariance;
        f.npar = npar;
        f
    }

    /// Barycentric J2000 states of `indices` at the TDB Julian dates `targets`, row-major
    /// `[object][target]`, integrated (N-body, or two-body with `opts.method`) from each
    /// object's epoch. NaN where an integration failed.
    pub fn states_at(&self, indices: &[usize], targets: &[f64], kernel: &SpiceKernel, opts: &BatchOptions) -> Vec<[f64; 6]> {
        self.integrate(indices, targets, kernel, opts, false)
    }

    /// [`Catalog::states_at`], starting each object from its snapshot state when that is closer
    /// in time to the targets than its own epoch.
    pub(crate) fn states_at_nearest(&self, indices: &[usize], targets: &[f64], kernel: &SpiceKernel, opts: &BatchOptions) -> Vec<[f64; 6]> {
        self.integrate(indices, targets, kernel, opts, true)
    }

    fn integrate(&self, indices: &[usize], targets: &[f64], kernel: &SpiceKernel, opts: &BatchOptions, nearest: bool) -> Vec<[f64; 6]> {
        let m = targets.len();
        let mut out = vec![[f64::NAN; 6]; indices.len() * m];
        if m == 0 || indices.is_empty() {
            return out;
        }
        let mid = 0.5 * (targets.iter().copied().fold(f64::INFINITY, f64::min) + targets.iter().copied().fold(f64::NEG_INFINITY, f64::max));
        let snap = nearest && self.has_snapshot();
        let start = |i: usize| -> (f64, [f64; 6]) {
            if snap && self.snapshot[i][0].is_finite() && (self.snapshot_epoch - mid).abs() < (self.orbits.epochs[i] - mid).abs() {
                (self.snapshot_epoch, self.snapshot[i])
            } else {
                (self.orbits.epochs[i], self.orbits.states[i])
            }
        };
        // The orbits with finite starting states, unnamed, carrying their non-gravitational
        // parameters.
        let mut ok = Vec::with_capacity(indices.len());
        let mut pop = Population::new(ReferencePlane::J2000, Origin::SSB);
        for (k, &i) in indices.iter().enumerate() {
            let (t0, s) = start(i);
            if !(t0.is_finite() && s.iter().all(|x| x.is_finite())) {
                continue;
            }
            let ng = self.nongrav[i];
            let properties = ng.iter().any(|&a| a != 0.0).then(|| Properties { nongrav: Some(ng), ..Default::default() });
            pop.states.push(s);
            pop.epochs.push(t0);
            pop.names.push(String::new());
            pop.properties.push(properties);
            ok.push(k);
        }
        match states_at(&pop, targets, kernel, opts) {
            Ok(s) => {
                for (j, &k) in ok.iter().enumerate() {
                    out[k * m..(k + 1) * m].copy_from_slice(&s[j * m..(j + 1) * m]);
                }
            }
            Err(_) if pop.len() > 1 => {
                // One bad orbit (e.g. out of the ephemeris' range) fails its whole chunk: retry
                // one by one.
                use rayon::prelude::*;
                let single = BatchOptions { parallel: false, ..opts.clone() };
                let each: Vec<Option<Vec<[f64; 6]>>> = (0..pop.len()).into_par_iter().map(|j| states_at(&pop.select(&[j]), targets, kernel, &single).ok()).collect();
                for (j, &k) in ok.iter().enumerate() {
                    if let Some(s) = &each[j] {
                        out[k * m..(k + 1) * m].copy_from_slice(s);
                    }
                }
            }
            Err(_) => {}
        }
        out
    }

    /// Integrate every orbit to `epoch` (TDB Julian date) and keep the states as the catalog's
    /// snapshot. N-body by default (`opts.method`); objects are processed in blocks of
    /// `block` to bound memory.
    pub fn make_snapshot(&mut self, epoch: f64, kernel: &SpiceKernel, opts: &BatchOptions, block: usize) {
        let n = self.len();
        let block = block.max(1);
        // Orbits already at the epoch are copied.
        let mut snap: Vec<[f64; 6]> = (0..n).map(|i| if self.orbits.epochs[i] == epoch { self.orbits.states[i] } else { [f64::NAN; 6] }).collect();
        let todo: Vec<usize> = (0..n).filter(|&i| self.orbits.epochs[i] != epoch).collect();
        for idx in todo.chunks(block) {
            for (k, s) in idx.iter().zip(self.states_at(idx, &[epoch], kernel, opts)) {
                snap[*k] = s;
            }
        }
        self.snapshot = snap;
        self.snapshot_epoch = epoch;
    }

    /// Write the catalog (with its snapshot) to a binary file.
    pub fn save(&self, path: &Path) -> Result<(), BoxError> {
        let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
        let snap = self.has_snapshot();
        w.write_all(b"SRCAT\x02")?;
        w.write_all(&(self.len() as u64).to_le_bytes())?;
        w.write_all(&self.snapshot_epoch.to_le_bytes())?;
        w.write_all(&[snap as u8])?;
        let f = |w: &mut BufWriter<File>, x: f64| w.write_all(&x.to_le_bytes());
        for i in 0..self.len() {
            let b = self.orbits.names[i].as_bytes();
            w.write_all(&(b.len() as u32).to_le_bytes())?;
            w.write_all(b)?;
            f(&mut w, self.orbits.epochs[i])?;
            for x in self.orbits.states[i] {
                f(&mut w, x)?;
            }
            for x in self.nongrav[i] {
                f(&mut w, x)?;
            }
            let bits = self.fit_nongrav[i].iter().enumerate().fold(0u8, |acc, (k, &b)| acc | ((b as u8) << k));
            w.write_all(&[bits])?;
            f(&mut w, self.h[i])?;
            f(&mut w, self.g[i])?;
            f(&mut w, self.u[i])?;
            f(&mut w, self.last_obs[i])?;
            let c = &self.covariance[i];
            w.write_all(&(c.len() as u32).to_le_bytes())?;
            for &x in c {
                f(&mut w, x)?;
            }
            if snap {
                for x in self.snapshot[i] {
                    f(&mut w, x)?;
                }
            }
        }
        w.flush()?;
        Ok(())
    }

    /// Read a catalog written by [`Catalog::save`].
    pub fn load(path: &Path) -> Result<Catalog, BoxError> {
        let mut r = BufReader::with_capacity(1 << 20, File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?);
        let mut magic = [0u8; 6];
        r.read_exact(&mut magic)?;
        if &magic != b"SRCAT\x02" {
            return Err(format!("{} is not a spacerocks catalog", path.display()).into());
        }
        let mut b8 = [0u8; 8];
        let mut b4 = [0u8; 4];
        let mut b1 = [0u8; 1];
        r.read_exact(&mut b8)?;
        let n = u64::from_le_bytes(b8) as usize;
        r.read_exact(&mut b8)?;
        let snapshot_epoch = f64::from_le_bytes(b8);
        r.read_exact(&mut b1)?;
        let snap = b1[0] != 0;
        let mut cat = Catalog { snapshot_epoch, ..Default::default() };
        let mut f = |r: &mut BufReader<File>| -> Result<f64, BoxError> {
            r.read_exact(&mut b8)?;
            Ok(f64::from_le_bytes(b8))
        };
        for _ in 0..n {
            r.read_exact(&mut b4)?;
            let mut name = vec![0u8; u32::from_le_bytes(b4) as usize];
            r.read_exact(&mut name)?;
            let epoch = f(&mut r)?;
            let mut state = [0.0; 6];
            for x in state.iter_mut() {
                *x = f(&mut r)?;
            }
            let mut nongrav = [0.0; 3];
            for x in nongrav.iter_mut() {
                *x = f(&mut r)?;
            }
            r.read_exact(&mut b1)?;
            let fit_nongrav = [b1[0] & 1 != 0, b1[0] & 2 != 0, b1[0] & 4 != 0];
            let (h, g, u, last_obs) = (f(&mut r)?, f(&mut r)?, f(&mut r)?, f(&mut r)?);
            r.read_exact(&mut b4)?;
            let nc = u32::from_le_bytes(b4) as usize;
            let mut cov = Vec::with_capacity(nc);
            for _ in 0..nc {
                cov.push(f(&mut r)?);
            }
            cat.push(String::from_utf8(name)?, epoch, state, nongrav, fit_nongrav, cov, h, g, u, last_obs);
            if snap {
                let mut s = [0.0; 6];
                for x in s.iter_mut() {
                    *x = f(&mut r)?;
                }
                cat.snapshot.push(s);
            }
        }
        if !snap {
            cat.snapshot_epoch = f64::NAN;
        }
        Ok(cat)
    }
}

