//! Observations for orbit fitting, stored as parallel arrays.

use std::collections::HashMap;

use nalgebra::Vector3;
use rayon::prelude::*;

use crate::observing::{Observation, Observatory};
use crate::spice::SpiceKernel;
use crate::time::Time;

/// One arcsecond in radians.
pub const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// Default astrometric uncertainty (radians): one arcsecond as layup rounds it, 1/206265.
pub const DEFAULT_SIGMA: f64 = 1.0 / 206265.0;

/// Default uncertainty of a sky-motion rate (radians/day): 24"/day (1"/hour), as layup.
pub const DEFAULT_RATE_SIGMA: f64 = 24.0 / 206265.0;

/// Default uncertainty of a radar delay (days): 1 µs, as layup.
pub const DEFAULT_DELAY_SIGMA: f64 = 1e-6 / 86400.0;

/// Observations of one object: one entry per detection in each array.
///
/// Every detection has an `epoch` (TDB Julian date; for radar, the receive time) and an
/// `observer` position (barycentric J2000, AU). What was measured is given by the other columns,
/// NaN where a detection has no such measurement, so optical astrometry, streaks and radar can
/// be mixed freely (as in an ADES table):
///
/// - `ra`, `dec` (radians, ICRF) with `sigma_ra`, `sigma_dec`. `sigma_ra` is the on-sky
///   uncertainty, of `ra * cos(dec)`, as ADES `rmsRA`.
/// - Optional sky-motion rates `ra_rate`, `dec_rate` (radians/day) with `sigma_ra_rate`,
///   `sigma_dec_rate`. `ra_rate` is the on-sky (great-circle) rate `cos(dec) dRA/dt`, as ADES
///   `raRate`.
/// - Optional radar: `delay` (round-trip light time, days) and `doppler` (round-trip range rate
///   d(delay)/dt times c, AU/day; a receding object has a positive value), with `sigma_delay`,
///   `sigma_doppler`.
///
/// Rates and radar need `observer_velocity` (AU/day). For radar, the transmitter is placed at the
/// transmit time (receive time minus the round-trip light time) by, in order of preference:
///
/// 1. `transmitter`: its barycentric state (position and velocity) at the transmit time, given
///    directly;
/// 2. `transmitter_site`: its Earth-fixed (ITRF93) position in AU, from which the state is
///    computed exactly at the model's transmit time. For a ground station this is the right
///    choice, and for monostatic radar it is the receiving station's own position;
/// 3. otherwise, the receiving station extrapolated back with `observer_velocity` and
///    `observer_acceleration` (AU/day²; zero if absent), as layup does. Over a round trip of a
///    few minutes this leaves ~0.1 m/s of Earth-rotation error in the transmitter's velocity
///    (~1 Hz of Doppler at S band), so prefer `transmitter_site`.
///
/// NaN entries in `transmitter` or `transmitter_site` fall through to the next option. Optional
/// columns are either empty or as long as `epoch`.
#[derive(Debug, Clone, Default)]
pub struct Astrometry {
    pub epoch: Vec<f64>,
    pub ra: Vec<f64>,
    pub dec: Vec<f64>,
    pub sigma_ra: Vec<f64>,
    pub sigma_dec: Vec<f64>,
    pub observer: Vec<[f64; 3]>,
    pub observer_velocity: Vec<[f64; 3]>,
    pub observer_acceleration: Vec<[f64; 3]>,
    pub ra_rate: Vec<f64>,
    pub dec_rate: Vec<f64>,
    pub sigma_ra_rate: Vec<f64>,
    pub sigma_dec_rate: Vec<f64>,
    pub delay: Vec<f64>,
    pub doppler: Vec<f64>,
    pub sigma_delay: Vec<f64>,
    pub sigma_doppler: Vec<f64>,
    pub transmitter: Vec<[f64; 6]>,
    pub transmitter_site: Vec<[f64; 3]>,
}

/// What a residual row measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowKind {
    Ra = 0,
    Dec = 1,
    RaRate = 2,
    DecRate = 3,
    Delay = 4,
    Doppler = 5,
}

fn positive_or(s: f64, default: f64) -> f64 {
    if s.is_finite() && s > 0.0 { s } else { default }
}

impl Astrometry {
    /// Optical astrometry with the default uncertainty ([`DEFAULT_SIGMA`], 1") on each axis.
    pub fn new(epoch: Vec<f64>, ra: Vec<f64>, dec: Vec<f64>, observer: Vec<[f64; 3]>) -> Result<Astrometry, String> {
        let n = epoch.len();
        let a = Astrometry { sigma_ra: vec![DEFAULT_SIGMA; n], sigma_dec: vec![DEFAULT_SIGMA; n], epoch, ra, dec, observer, ..Default::default() };
        a.check()?;
        Ok(a)
    }

    /// Replace the astrometric uncertainties (radians). Non-finite or non-positive entries get
    /// [`DEFAULT_SIGMA`].
    pub fn with_sigma(mut self, sigma_ra: &[f64], sigma_dec: &[f64]) -> Result<Astrometry, String> {
        self.sigma_ra = sigma_ra.iter().map(|&s| positive_or(s, DEFAULT_SIGMA)).collect();
        self.sigma_dec = sigma_dec.iter().map(|&s| positive_or(s, DEFAULT_SIGMA)).collect();
        self.check()?;
        Ok(self)
    }

    /// Add observer velocities (AU/day), needed for rates and radar.
    pub fn with_observer_velocity(mut self, velocity: Vec<[f64; 3]>) -> Result<Astrometry, String> {
        self.observer_velocity = velocity;
        self.check()?;
        Ok(self)
    }

    /// Add observer accelerations (AU/day²), used by monostatic radar.
    pub fn with_observer_acceleration(mut self, acceleration: Vec<[f64; 3]>) -> Result<Astrometry, String> {
        self.observer_acceleration = acceleration;
        self.check()?;
        Ok(self)
    }

    /// Add sky-motion rates (radians/day; `ra_rate` on-sky, i.e. times cos Dec). NaN rates mark
    /// detections without them; non-finite or non-positive uncertainties get
    /// [`DEFAULT_RATE_SIGMA`].
    pub fn with_rates(mut self, ra_rate: Vec<f64>, dec_rate: Vec<f64>, sigma_ra_rate: &[f64], sigma_dec_rate: &[f64]) -> Result<Astrometry, String> {
        self.ra_rate = ra_rate;
        self.dec_rate = dec_rate;
        self.sigma_ra_rate = sigma_ra_rate.iter().map(|&s| positive_or(s, DEFAULT_RATE_SIGMA)).collect();
        self.sigma_dec_rate = sigma_dec_rate.iter().map(|&s| positive_or(s, DEFAULT_RATE_SIGMA)).collect();
        self.check()?;
        Ok(self)
    }

    /// Add radar delay (days) and Doppler (round-trip range rate, AU/day) measurements; NaN
    /// marks detections without one. A missing delay uncertainty gets [`DEFAULT_DELAY_SIGMA`];
    /// every Doppler measurement needs an uncertainty.
    pub fn with_radar(mut self, delay: Vec<f64>, doppler: Vec<f64>, sigma_delay: &[f64], sigma_doppler: &[f64]) -> Result<Astrometry, String> {
        self.sigma_delay = sigma_delay.iter().map(|&s| positive_or(s, DEFAULT_DELAY_SIGMA)).collect();
        self.sigma_doppler = sigma_doppler.to_vec();
        if doppler.iter().zip(sigma_doppler).any(|(d, s)| d.is_finite() && !(s.is_finite() && *s > 0.0)) {
            return Err("every Doppler measurement needs a positive uncertainty".into());
        }
        self.delay = delay;
        self.doppler = doppler;
        self.check()?;
        Ok(self)
    }

    /// Add the transmitting antenna's barycentric state (AU, AU/day) at the transmit time, NaN
    /// for detections whose receiving station transmitted.
    pub fn with_transmitter(mut self, transmitter: Vec<[f64; 6]>) -> Result<Astrometry, String> {
        self.transmitter = transmitter;
        self.check()?;
        Ok(self)
    }

    /// Add the Earth-fixed (ITRF93) position of the transmitting antenna (AU), NaN where
    /// unknown. See [`crate::observing::Observatory::earth_fixed_position`].
    pub fn with_transmitter_site(mut self, site: Vec<[f64; 3]>) -> Result<Astrometry, String> {
        self.transmitter_site = site;
        self.check()?;
        Ok(self)
    }

    /// Astrometry from [`Observation`]s (astrometry or streaks). The uncertainties come from
    /// the observations' covariances when present.
    pub fn from_observations(observations: &[Observation]) -> Result<Astrometry, Box<dyn std::error::Error>> {
        let mut a = Astrometry::default();
        let rates = observations.iter().any(|o| o.ra_rate().is_some());
        for o in observations {
            let mut obs = o.observer.clone();
            if obs.origin.name().to_uppercase() != "SSB" {
                return Err(format!("observer positions must be barycentric (got origin {})", obs.origin.name()).into());
            }
            if obs.reference_plane.as_str() != "J2000" {
                obs.change_reference_plane("J2000")?;
            }
            let dec = o.dec();
            a.epoch.push(o.epoch.tdb().jd());
            a.ra.push(o.ra());
            a.dec.push(dec);
            let cov = o.covariance();
            let sd = |k: usize| cov.as_ref().filter(|c| c.nrows() > k).map(|c| c[(k, k)].sqrt()).unwrap_or(f64::NAN);
            a.sigma_ra.push(positive_or(sd(0) * dec.cos(), DEFAULT_SIGMA));
            a.sigma_dec.push(positive_or(sd(1), DEFAULT_SIGMA));
            a.observer.push([obs.position.x, obs.position.y, obs.position.z]);
            if rates {
                let v = obs.velocity.ok_or("rates need the observer's velocity")?;
                a.observer_velocity.push([v.x, v.y, v.z]);
                // Observation rates are coordinate rates (dRA/dt); the fit uses on-sky rates.
                a.ra_rate.push(o.ra_rate().map(|r| r * dec.cos()).unwrap_or(f64::NAN));
                a.dec_rate.push(o.dec_rate().unwrap_or(f64::NAN));
                a.sigma_ra_rate.push(positive_or(sd(2) * dec.cos(), DEFAULT_RATE_SIGMA));
                a.sigma_dec_rate.push(positive_or(sd(3), DEFAULT_RATE_SIGMA));
            }
        }
        a.check()?;
        Ok(a)
    }

    fn check(&self) -> Result<(), String> {
        let n = self.epoch.len();
        let required = [("ra", self.ra.len()), ("dec", self.dec.len()), ("sigma_ra", self.sigma_ra.len()), ("sigma_dec", self.sigma_dec.len()), ("observer", self.observer.len())];
        for (name, m) in required {
            if m != n {
                return Err(format!("{} has {} entries for {} epochs", name, m, n));
            }
        }
        let optional = [
            ("observer_velocity", self.observer_velocity.len()),
            ("observer_acceleration", self.observer_acceleration.len()),
            ("ra_rate", self.ra_rate.len()),
            ("dec_rate", self.dec_rate.len()),
            ("sigma_ra_rate", self.sigma_ra_rate.len()),
            ("sigma_dec_rate", self.sigma_dec_rate.len()),
            ("delay", self.delay.len()),
            ("doppler", self.doppler.len()),
            ("sigma_delay", self.sigma_delay.len()),
            ("sigma_doppler", self.sigma_doppler.len()),
            ("transmitter", self.transmitter.len()),
            ("transmitter_site", self.transmitter_site.len()),
        ];
        for (name, m) in optional {
            if m != 0 && m != n {
                return Err(format!("{} has {} entries for {} epochs", name, m, n));
            }
        }
        if (!self.ra_rate.is_empty() || !self.delay.is_empty()) && self.observer_velocity.is_empty() {
            return Err("rates and radar need observer velocities".into());
        }
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.epoch.len()
    }

    pub fn is_empty(&self) -> bool {
        self.epoch.is_empty()
    }

    /// Whether detection `i` has optical astrometry (finite RA and Dec).
    #[inline]
    pub fn has_optical(&self, i: usize) -> bool {
        self.ra[i].is_finite() && self.dec[i].is_finite()
    }

    #[inline]
    pub fn has_rates(&self, i: usize) -> bool {
        !self.ra_rate.is_empty() && self.ra_rate[i].is_finite() && self.dec_rate[i].is_finite() && self.has_optical(i)
    }

    #[inline]
    pub fn has_delay(&self, i: usize) -> bool {
        !self.delay.is_empty() && self.delay[i].is_finite()
    }

    #[inline]
    pub fn has_doppler(&self, i: usize) -> bool {
        !self.doppler.is_empty() && self.doppler[i].is_finite()
    }

    /// The residual rows detection `i` contributes, in order.
    pub fn rows(&self, i: usize) -> impl Iterator<Item = RowKind> {
        let opt = self.has_optical(i);
        let rates = self.has_rates(i);
        let (delay, doppler) = (self.has_delay(i), self.has_doppler(i));
        [(opt, RowKind::Ra), (opt, RowKind::Dec), (rates, RowKind::RaRate), (rates, RowKind::DecRate), (delay, RowKind::Delay), (doppler, RowKind::Doppler)]
            .into_iter()
            .filter(|(on, _)| *on)
            .map(|(_, k)| k)
    }

    /// Total number of residual rows.
    pub fn n_rows(&self) -> usize {
        (0..self.len()).map(|i| self.rows(i).count()).sum()
    }

    /// Uncertainty of the measurement behind a residual row.
    #[inline]
    pub fn sigma(&self, i: usize, kind: RowKind) -> f64 {
        match kind {
            RowKind::Ra => self.sigma_ra[i],
            RowKind::Dec => self.sigma_dec[i],
            RowKind::RaRate => self.sigma_ra_rate[i],
            RowKind::DecRate => self.sigma_dec_rate[i],
            RowKind::Delay => self.sigma_delay[i],
            RowKind::Doppler => self.sigma_doppler[i],
        }
    }

    /// Indices of the detections with optical astrometry.
    pub fn optical(&self) -> Vec<usize> {
        (0..self.len()).filter(|&i| self.has_optical(i)).collect()
    }

    /// The detections at `idx`, in that order.
    pub fn subset(&self, idx: &[usize]) -> Astrometry {
        fn pick<T: Copy>(v: &[T], idx: &[usize]) -> Vec<T> {
            if v.is_empty() { Vec::new() } else { idx.iter().map(|&i| v[i]).collect() }
        }
        Astrometry {
            epoch: pick(&self.epoch, idx),
            ra: pick(&self.ra, idx),
            dec: pick(&self.dec, idx),
            sigma_ra: pick(&self.sigma_ra, idx),
            sigma_dec: pick(&self.sigma_dec, idx),
            observer: pick(&self.observer, idx),
            observer_velocity: pick(&self.observer_velocity, idx),
            observer_acceleration: pick(&self.observer_acceleration, idx),
            ra_rate: pick(&self.ra_rate, idx),
            dec_rate: pick(&self.dec_rate, idx),
            sigma_ra_rate: pick(&self.sigma_ra_rate, idx),
            sigma_dec_rate: pick(&self.sigma_dec_rate, idx),
            delay: pick(&self.delay, idx),
            doppler: pick(&self.doppler, idx),
            sigma_delay: pick(&self.sigma_delay, idx),
            sigma_doppler: pick(&self.sigma_doppler, idx),
            transmitter: pick(&self.transmitter, idx),
            transmitter_site: pick(&self.transmitter_site, idx),
        }
    }

    /// A 64-bit key for detection `i`: an FNV-1a hash of every value the fit uses (epoch,
    /// measurements, uncertainties, observer and transmitter states). Equal detections have equal
    /// keys, on any machine. Used to tell which detections are new since a previous fit.
    pub fn detection_key(&self, i: usize) -> u64 {
        const PRIME: u64 = 0x100000001b3;
        let mut h: u64 = 0xcbf29ce484222325;
        let mut mix = |v: f64| {
            // one bit pattern for every NaN, and for +0 and -0
            let b = if v.is_nan() { f64::NAN.to_bits() } else if v == 0.0 { 0 } else { v.to_bits() };
            for byte in b.to_le_bytes() {
                h ^= byte as u64;
                h = h.wrapping_mul(PRIME);
            }
        };
        let get = |v: &Vec<f64>| v.get(i).copied().unwrap_or(f64::NAN);
        for v in [&self.epoch, &self.ra, &self.dec, &self.sigma_ra, &self.sigma_dec, &self.ra_rate, &self.dec_rate,
                  &self.sigma_ra_rate, &self.sigma_dec_rate, &self.delay, &self.doppler, &self.sigma_delay, &self.sigma_doppler] {
            mix(get(v));
        }
        for v in [&self.observer, &self.observer_velocity, &self.observer_acceleration, &self.transmitter_site] {
            for c in v.get(i).copied().unwrap_or([f64::NAN; 3]) {
                mix(c);
            }
        }
        for c in self.transmitter.get(i).copied().unwrap_or([f64::NAN; 6]) {
            mix(c);
        }
        h
    }

    /// Keys of all detections ([`Astrometry::detection_key`]).
    pub fn detection_keys(&self) -> Vec<u64> {
        (0..self.len()).map(|i| self.detection_key(i)).collect()
    }

    /// A fingerprint of the whole set of detections, independent of their order (layup's
    /// `obs_hash`, over the fit's inputs rather than the reported columns).
    pub fn fingerprint(&self) -> u64 {
        fingerprint(&self.detection_keys())
    }

    /// Indices that sort the detections by epoch (stable).
    pub fn time_order(&self) -> Vec<usize> {
        let mut idx: Vec<usize> = (0..self.len()).collect();
        idx.sort_by(|&a, &b| self.epoch[a].total_cmp(&self.epoch[b]));
        idx
    }

    /// Unit vector towards the detection.
    #[inline]
    pub fn rho_hat(&self, i: usize) -> Vector3<f64> {
        let (sd, cd) = self.dec[i].sin_cos();
        let (sa, ca) = self.ra[i].sin_cos();
        Vector3::new(cd * ca, cd * sa, sd)
    }

    /// Unit vectors of the tangent plane at the detection: along increasing RA (`a`) and Dec
    /// (`d`), built from the unit vector as in layup.
    #[inline]
    pub fn tangent_basis(&self, i: usize) -> (Vector3<f64>, Vector3<f64>) {
        let r = self.rho_hat(i);
        let cd = (r.x * r.x + r.y * r.y).sqrt();
        (Vector3::new(-r.y / cd, r.x / cd, 0.0), Vector3::new(-r.z * r.x / cd, -r.z * r.y / cd, cd))
    }

    #[inline]
    pub fn observer(&self, i: usize) -> Vector3<f64> {
        Vector3::from(self.observer[i])
    }

    #[inline]
    pub fn observer_velocity(&self, i: usize) -> Vector3<f64> {
        self.observer_velocity.get(i).map(|v| Vector3::from(*v)).unwrap_or_else(Vector3::zeros)
    }

    #[inline]
    pub fn observer_acceleration(&self, i: usize) -> Vector3<f64> {
        self.observer_acceleration.get(i).map(|v| Vector3::from(*v)).unwrap_or_else(Vector3::zeros)
    }

    /// Earth-fixed position of the transmitting antenna for radar detection `i`, if given.
    #[inline]
    pub fn transmitter_site(&self, i: usize) -> Option<[f64; 3]> {
        self.transmitter_site.get(i).copied().filter(|p| p[0].is_finite())
    }

    /// Transmitter position and velocity for radar detection `i`, if a separate one was given.
    #[inline]
    pub fn transmitter(&self, i: usize) -> Option<(Vector3<f64>, Vector3<f64>)> {
        let t = self.transmitter.get(i)?;
        if t[0].is_finite() {
            Some((Vector3::new(t[0], t[1], t[2]), Vector3::new(t[3], t[4], t[5])))
        } else {
            None
        }
    }
}

fn sites<'a>(codes: &'a [String]) -> Result<HashMap<&'a str, Observatory>, String> {
    let mut sites = HashMap::new();
    for c in codes {
        if !sites.contains_key(c.as_str()) {
            let o = Observatory::from_obscode(c).map_err(|e| format!("observatory '{}': {}", c, e))?;
            sites.insert(c.as_str(), o);
        }
    }
    Ok(sites)
}

/// Barycentric J2000 states (position AU, velocity AU/day) of the observatories with MPC codes
/// `codes` at the TDB Julian dates `epochs`.
pub fn observer_states(codes: &[String], epochs: &[f64], kernel: &SpiceKernel) -> Result<Vec<[f64; 6]>, String> {
    if codes.len() != epochs.len() {
        return Err(format!("{} observatory codes for {} epochs", codes.len(), epochs.len()));
    }
    let sites = sites(codes)?;
    codes
        .par_iter()
        .zip(epochs.par_iter())
        .map(|(c, &jd)| {
            let t = Time::new(jd, "tdb", "jd").map_err(|e| e.to_string())?;
            let o = sites[c.as_str()].at(&t, "J2000", "SSB", kernel).map_err(|e| format!("observatory '{}': {}", c, e))?;
            let v = o.velocity.unwrap_or_else(Vector3::zeros);
            Ok([o.position.x, o.position.y, o.position.z, v.x, v.y, v.z])
        })
        .collect()
}

/// Barycentric J2000 positions (AU) of the observatories with MPC codes `codes` at the TDB Julian
/// dates `epochs`.
pub fn observer_positions(codes: &[String], epochs: &[f64], kernel: &SpiceKernel) -> Result<Vec<[f64; 3]>, String> {
    Ok(observer_states(codes, epochs, kernel)?.into_iter().map(|s| [s[0], s[1], s[2]]).collect())
}

/// Barycentric accelerations (AU/day²) of the observatories, by central differences of their
/// velocities over ±2 s (as layup).
pub fn observer_accelerations(codes: &[String], epochs: &[f64], kernel: &SpiceKernel) -> Result<Vec<[f64; 3]>, String> {
    let dt = 2.0 / 86400.0;
    let plus: Vec<f64> = epochs.iter().map(|t| t + dt).collect();
    let minus: Vec<f64> = epochs.iter().map(|t| t - dt).collect();
    let (p, m) = (observer_states(codes, &plus, kernel)?, observer_states(codes, &minus, kernel)?);
    Ok(p.iter().zip(&m).map(|(p, m)| [(p[3] - m[3]) / (2.0 * dt), (p[4] - m[4]) / (2.0 * dt), (p[5] - m[5]) / (2.0 * dt)]).collect())
}

/// WGS84 ellipsoid: semi-major axis (m) and flattening.
const WGS84_A_M: f64 = 6378137.0;
const WGS84_F: f64 = 1.0 / 298.257223563;

/// Earth-fixed (ITRF93) position in AU of a point at WGS84 geodetic east longitude and latitude
/// (degrees) and height above the ellipsoid (m), as layup places a roving observer.
pub fn geodetic_to_earth_fixed(lon_deg: f64, lat_deg: f64, height_m: f64) -> [f64; 3] {
    let e2 = WGS84_F * (2.0 - WGS84_F);
    let (lon, lat) = (lon_deg * std::f64::consts::PI / 180.0, lat_deg * std::f64::consts::PI / 180.0);
    let (sl, cl) = lat.sin_cos();
    let n = WGS84_A_M / (1.0 - e2 * sl * sl).sqrt();
    let rho_cos = (n + height_m) * cl * crate::data::constants::M_TO_AU;
    let rho_sin = (n * (1.0 - e2) + height_m) * sl * crate::data::constants::M_TO_AU;
    [rho_cos * lon.cos(), rho_cos * lon.sin(), rho_sin]
}

/// Barycentric J2000 state (AU, AU/day) at TDB Julian date `jd_tdb` of an observer whose position
/// is given with the detection, in ADES form (layup's `populate_observatory`): `sys` is
/// `"ICRF_KM"` or `"ICRF_AU"` (a geocentric ICRF position in km or AU, with an optional
/// geocentric velocity in km/s or AU/day; without one the observer moves with the Earth's
/// center), or `"WGS84"` (`pos` = east longitude and geodetic latitude in degrees and height in
/// m, e.g. a roving observer, which then turns with the Earth). The center `ctr` must be 399
/// (the Earth), as layup requires.
pub fn ades_observer_state(sys: &str, ctr: i64, pos: [f64; 3], vel: Option<[f64; 3]>, jd_tdb: f64, kernel: &SpiceKernel) -> Result<[f64; 6], String> {
    if ctr != 399 {
        return Err(format!("observer center {} is not supported (use 399, the Earth)", ctr));
    }
    if !pos.iter().all(|v| v.is_finite()) {
        return Err("observer position is not finite".into());
    }
    let km = crate::data::constants::KM_TO_AU;
    let (scale_p, scale_v) = match sys {
        "WGS84" => {
            let ef = geodetic_to_earth_fixed(pos[0], pos[1], pos[2]);
            return crate::observing::earth_fixed_state(&ef, jd_tdb, kernel).map_err(|e| e.to_string());
        }
        "ICRF_KM" => (km, km * 86400.0),
        "ICRF_AU" => (1.0, 1.0),
        other => return Err(format!("observer frame '{}' is not supported (use ICRF_KM, ICRF_AU or WGS84)", other)),
    };
    let e = kernel.state_au(399, 0, jd_tdb).map_err(|e| e.to_string())?;
    let v = vel.filter(|v| v.iter().all(|x| x.is_finite())).unwrap_or([0.0; 3]);
    Ok([
        e[0] + pos[0] * scale_p,
        e[1] + pos[1] * scale_p,
        e[2] + pos[2] * scale_p,
        e[3] + v[0] * scale_v,
        e[4] + v[1] * scale_v,
        e[5] + v[2] * scale_v,
    ])
}

/// Position of an occulting object from ADES occultation astrometry: the occulted star's RA and
/// Dec (`raStar`, `decStar`) and the object's offset from it (`deltaRA`, which includes cos Dec,
/// and `deltaDec`), all in radians. The offset is applied on the tangent plane at the star (the
/// exact form of layup's intended `raStar + deltaRA / cos(decStar)`, `decStar + deltaDec`; they
/// differ at second order in the offset, ~offset² tan(Dec) / 2: under 0.3 mas for offsets of 2").
pub fn occultation_radec(ra_star: f64, dec_star: f64, delta_ra: f64, delta_dec: f64) -> (f64, f64) {
    let (sd, cd) = dec_star.sin_cos();
    let den = cd - delta_dec * sd;
    let ra = (ra_star + delta_ra.atan2(den)).rem_euclid(std::f64::consts::TAU);
    let dec = (sd + delta_dec * cd).atan2((delta_ra * delta_ra + den * den).sqrt());
    (ra, dec)
}

/// Order-independent fingerprint of a set of detection keys: FNV-1a over the sorted keys.
pub fn fingerprint(keys: &[u64]) -> u64 {
    let mut k = keys.to_vec();
    k.sort_unstable();
    let mut h: u64 = 0xcbf29ce484222325;
    for key in k {
        for byte in key.to_le_bytes() {
            h ^= byte as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}
