//! Reading the MPC's orbit catalogs.
//!
//! - `mpcorb_extended.json(.gz)` and `MPCORB.DAT(.gz)`: osculating heliocentric ecliptic
//!   elements (J2000) of every minor planet, at an epoch in TT. No covariance, but each orbit has
//!   the MPC's uncertainty parameter U (0–9).
//! - `mpc_orb` JSON (the MPC's `get-orb` API): one object's orbit, with the covariance of its
//!   heliocentric ecliptic Cartesian state ("CAR").

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use nalgebra::{Matrix3, Vector3};
use serde::Deserialize;

use crate::constants::ROTATION_ECLIPJ2000;
use crate::time::conversions::tt_to_tdb;
use crate::transforms::calc_true_anomaly_from_mean_anomaly;
use crate::Origin;

type BoxError = Box<dyn std::error::Error + Send + Sync>;

/// URL of the MPC's extended JSON catalog of all minor-planet orbits.
pub const MPCORB_JSON_URL: &str = "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz";
/// URL of the MPC's fixed-width catalog of all minor-planet orbits.
pub const MPCORB_DAT_URL: &str = "https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz";

/// One orbit read from an MPC catalog: heliocentric, ecliptic J2000 elements at `epoch_tt`.
#[derive(Debug, Clone)]
pub struct MpcElements {
    pub name: String,
    /// TT Julian date of the elements.
    pub epoch_tt: f64,
    /// Semimajor axis (AU), eccentricity, and inclination, ascending node, argument of
    /// perihelion and mean anomaly (degrees).
    pub a: f64,
    pub e: f64,
    pub inc: f64,
    pub node: f64,
    pub peri: f64,
    pub m: f64,
    pub h: f64,
    pub g: f64,
    /// The MPC uncertainty parameter U (0–9); NaN if the catalog gives a letter or nothing.
    pub u: f64,
    /// Julian date of the last observation used in the orbit (NaN if not given).
    pub last_obs: f64,
}

/// The directory MPC catalogs are kept in: `$SPACEROCKS_MPC_DIR`, else `~/.spacerocks/mpc`.
pub fn mpc_dir() -> PathBuf {
    if let Ok(d) = std::env::var("SPACEROCKS_MPC_DIR") {
        return PathBuf::from(d);
    }
    dirs::home_dir().unwrap_or_default().join(".spacerocks").join("mpc")
}

/// Download `url` to `path` (via a temporary file, so an interrupted download leaves nothing
/// behind).
pub fn download(url: &str, path: &Path) -> Result<(), BoxError> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let client = reqwest::blocking::Client::builder().timeout(std::time::Duration::from_secs(3600)).build()?;
    let mut resp = client.get(url).send()?.error_for_status()?;
    let tmp = path.with_extension("part");
    {
        let mut f = File::create(&tmp)?;
        resp.copy_to(&mut f)?;
    }
    std::fs::rename(&tmp, path)?;
    Ok(())
}

/// The MPC's extended JSON catalog, downloaded to [`mpc_dir`] if it is not there yet (or if
/// `update`). Returns its path.
pub fn mpcorb_path(download_missing: bool, update: bool) -> Result<PathBuf, BoxError> {
    let path = mpc_dir().join("mpcorb_extended.json.gz");
    if update || !path.exists() {
        if !download_missing && !update {
            return Err(format!("{} does not exist (download it with download=True, or from {})", path.display(), MPCORB_JSON_URL).into());
        }
        download(MPCORB_JSON_URL, &path)?;
    }
    Ok(path)
}

fn open(path: &Path) -> Result<Box<dyn Read>, BoxError> {
    let f = File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let gz = path.extension().map(|x| x == "gz").unwrap_or(false);
    Ok(if gz { Box::new(GzDecoder::new(BufReader::with_capacity(1 << 20, f))) } else { Box::new(BufReader::with_capacity(1 << 20, f)) })
}

#[derive(Deserialize)]
#[allow(non_snake_case)]
struct JsonEntry {
    #[serde(default)]
    Number: Option<String>,
    #[serde(default)]
    Principal_desig: Option<String>,
    #[serde(default)]
    Name: Option<String>,
    Epoch: f64,
    a: f64,
    e: f64,
    i: f64,
    Node: f64,
    Peri: f64,
    M: f64,
    #[serde(default)]
    H: Option<f64>,
    #[serde(default)]
    G: Option<f64>,
    #[serde(default)]
    U: Option<serde_json::Value>,
    #[serde(default)]
    Last_obs: Option<String>,
}

/// "YYYY-MM-DD" or "YYYYMMDD" as a Julian date (0h).
fn date_jd(s: &str) -> f64 {
    let d: String = s.chars().filter(|c| c.is_ascii_digit()).collect();
    if d.len() != 8 {
        return f64::NAN;
    }
    match (d[0..4].parse::<i32>(), d[4..6].parse::<i32>(), d[6..8].parse::<f64>()) {
        (Ok(y), Ok(m), Ok(day)) => calendar_to_jd(y, m, day),
        _ => f64::NAN,
    }
}

fn parse_u(s: &str) -> f64 {
    match s.trim().parse::<u8>() {
        Ok(u) if u <= 9 => u as f64,
        _ => f64::NAN,
    }
}

/// Read `mpcorb_extended.json` (optionally gzipped). Numbered objects are named by their
/// number (e.g. "3666"), the others by their principal designation.
pub fn read_mpcorb_json(path: &Path) -> Result<Vec<MpcElements>, BoxError> {
    let entries: Vec<JsonEntry> = serde_json::from_reader(open(path)?)?;
    Ok(entries
        .into_iter()
        .map(|x| {
            let name = match (&x.Number, &x.Principal_desig, &x.Name) {
                (Some(n), _, _) => n.trim_matches(|c| c == '(' || c == ')').to_string(),
                (None, Some(d), _) => d.clone(),
                (None, None, Some(n)) => n.clone(),
                _ => String::new(),
            };
            let u = match &x.U {
                Some(serde_json::Value::String(s)) => parse_u(s),
                Some(serde_json::Value::Number(n)) => n.as_f64().filter(|u| (0.0..=9.0).contains(u)).unwrap_or(f64::NAN),
                _ => f64::NAN,
            };
            MpcElements {
                name,
                epoch_tt: x.Epoch,
                a: x.a,
                e: x.e,
                inc: x.i,
                node: x.Node,
                peri: x.Peri,
                m: x.M,
                h: x.H.unwrap_or(f64::NAN),
                g: x.G.unwrap_or(0.15),
                u,
                last_obs: x.Last_obs.as_deref().map(date_jd).unwrap_or(f64::NAN),
            }
        })
        .collect())
}

fn packed_digit(c: u8) -> Option<u32> {
    match c {
        b'0'..=b'9' => Some((c - b'0') as u32),
        b'A'..=b'Z' => Some((c - b'A') as u32 + 10),
        b'a'..=b'z' => Some((c - b'a') as u32 + 36),
        _ => None,
    }
}

/// A packed MPC date ("K2555" = 2025 May 5.0) as a Julian date.
pub fn unpack_epoch(s: &str) -> Option<f64> {
    let b = s.trim().as_bytes();
    if b.len() != 5 {
        return None;
    }
    let century = match b[0] {
        b'I' => 1800,
        b'J' => 1900,
        b'K' => 2000,
        b'L' => 2100,
        _ => return None,
    };
    let year = century + std::str::from_utf8(&b[1..3]).ok()?.parse::<i32>().ok()?;
    let month = packed_digit(b[3])? as i32;
    let day = packed_digit(b[4])? as f64;
    Some(calendar_to_jd(year, month, day))
}

/// Julian date of a Gregorian calendar date (day may be fractional).
fn calendar_to_jd(year: i32, month: i32, day: f64) -> f64 {
    let (y, m) = if month <= 2 { (year - 1, month + 12) } else { (year, month) };
    let a = y.div_euclid(100);
    let b = 2 - a + a.div_euclid(4);
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (m as f64 + 1.0)).floor() + day + b as f64 - 1524.5
}

/// Unpack an MPC packed designation: numbers ("03666" → "3666", "A0345" → "100345",
/// "~0000" → "620000") and provisional designations ("K07Tf8A" → "2007 TA418"). Anything else
/// is returned trimmed.
pub fn unpack_designation(s: &str) -> String {
    let t = s.trim();
    let b = t.as_bytes();
    if b.len() == 5 {
        if b[0] == b'~' {
            let mut n: u64 = 0;
            for &c in &b[1..] {
                let Some(d) = packed_digit(c) else { return t.to_string() };
                n = n * 62 + d as u64;
            }
            return (620_000 + n).to_string();
        }
        if let (Some(d), Ok(rest)) = (packed_digit(b[0]), t[1..].parse::<u64>()) {
            return (d as u64 * 10_000 + rest).to_string();
        }
    }
    if b.len() == 7 {
        let century = match b[0] {
            b'I' => 18,
            b'J' => 19,
            b'K' => 20,
            _ => return t.to_string(),
        };
        let cycle = packed_digit(b[4]).map(|d| d * 10).unwrap_or(0) + (b[5] as char).to_digit(10).unwrap_or(0);
        let num = if cycle == 0 { String::new() } else { cycle.to_string() };
        return format!("{}{} {}{}{}", century, &t[1..3], b[3] as char, b[6] as char, num);
    }
    t.to_string()
}

/// Read `MPCORB.DAT` (optionally gzipped): the header is skipped, and so are lines that do not
/// parse as orbits.
pub fn read_mpcorb_dat(path: &Path) -> Result<Vec<MpcElements>, BoxError> {
    let reader = BufReader::with_capacity(1 << 20, open(path)?);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.len() < 103 || !line.is_char_boundary(103) {
            continue;
        }
        let f = |a: usize, b: usize| line.get(a..b.min(line.len())).map(str::trim).unwrap_or("");
        let num = |a: usize, b: usize| f(a, b).parse::<f64>().ok();
        let (Some(epoch), Some(m), Some(peri), Some(node), Some(inc), Some(e), Some(a)) =
            (unpack_epoch(f(20, 25)), num(26, 35), num(37, 46), num(48, 57), num(59, 68), num(70, 79), num(92, 103))
        else {
            continue;
        };
        let name = unpack_designation(f(0, 7));
        out.push(MpcElements {
            name,
            epoch_tt: epoch,
            a,
            e,
            inc,
            node,
            peri,
            m,
            h: num(8, 13).unwrap_or(f64::NAN),
            g: num(14, 19).unwrap_or(0.15),
            u: parse_u(f(105, 106)),
            last_obs: date_jd(f(194, 202)),
        });
    }
    Ok(out)
}

/// Read an MPC catalog, by extension: `.json`/`.json.gz` or anything else as `MPCORB.DAT`.
pub fn read_mpcorb(path: &Path) -> Result<Vec<MpcElements>, BoxError> {
    let name = path.file_name().map(|s| s.to_string_lossy().to_lowercase()).unwrap_or_default();
    if name.ends_with(".json") || name.ends_with(".json.gz") {
        read_mpcorb_json(path)
    } else {
        read_mpcorb_dat(path)
    }
}

/// Matrix from the ecliptic J2000 frame to the J2000 equator, as `SpaceRock` uses it.
pub(crate) fn ecliptic_to_j2000() -> Matrix3<f64> {
    ROTATION_ECLIPJ2000.try_inverse().unwrap_or_else(Matrix3::identity)
}

/// Heliocentric state (AU, AU/day) in the elements' own frame from Keplerian elements (angles
/// in radians), for the Sun's GM. Works for any conic (a < 0 for hyperbolae).
pub fn elements_to_state(a: f64, e: f64, inc: f64, node: f64, peri: f64, mean_anomaly: f64) -> Option<(Vector3<f64>, Vector3<f64>)> {
    let mu = Origin::SUN.mu();
    let f = calc_true_anomaly_from_mean_anomaly(e, mean_anomaly).ok()?;
    let p = a * (1.0 - e * e);
    if !(p > 0.0) || !f.is_finite() {
        return None;
    }
    let r = p / (1.0 + e * f.cos());
    let pos_pf = Vector3::new(r * f.cos(), r * f.sin(), 0.0);
    let k = (mu / p).sqrt();
    let vel_pf = Vector3::new(-k * f.sin(), k * (e + f.cos()), 0.0);
    let (so, co) = node.sin_cos();
    let (si, ci) = inc.sin_cos();
    let (sw, cw) = peri.sin_cos();
    let rot = Matrix3::new(
        co * cw - so * sw * ci, -co * sw - so * cw * ci, so * si,
        so * cw + co * sw * ci, -so * sw + co * cw * ci, -co * si,
        sw * si, cw * si, ci,
    );
    Some((rot * pos_pf, rot * vel_pf))
}

impl MpcElements {
    /// TDB Julian date of the elements.
    pub fn epoch_tdb(&self) -> f64 {
        tt_to_tdb(self.epoch_tt)
    }

    /// Heliocentric J2000 (equatorial) state at [`MpcElements::epoch_tdb`].
    pub fn helio_j2000(&self) -> Option<[f64; 6]> {
        let d = std::f64::consts::PI / 180.0;
        let (p, v) = elements_to_state(self.a, self.e, self.inc * d, self.node * d, self.peri * d, self.m * d)?;
        let m = ecliptic_to_j2000();
        let (p, v) = (m * p, m * v);
        Some([p.x, p.y, p.z, v.x, v.y, v.z])
    }
}

/// An orbit from the MPC's `mpc_orb` JSON (e.g. from the `get-orb` API): heliocentric J2000
/// state with its 6x6 covariance.
#[derive(Debug, Clone)]
pub struct MpcOrb {
    pub name: String,
    pub epoch_tdb: f64,
    /// Heliocentric J2000 (equatorial) state.
    pub helio_j2000: [f64; 6],
    /// Row-major 6x6 covariance of `helio_j2000`.
    pub covariance: [f64; 36],
    pub h: f64,
    pub g: f64,
    pub u: f64,
    /// Non-gravitational parameters A1–A3 if the MPC fitted them (the covariance then does not
    /// include them).
    pub nongrav: [f64; 3],
}

/// Parse one `mpc_orb` JSON object (or a list with one, as the API wraps it).
pub fn parse_mpc_orb(v: &serde_json::Value) -> Result<MpcOrb, BoxError> {
    let mut v = v;
    while let Some(first) = v.as_array().and_then(|a| a.first()) {
        v = first;
    }
    if let Some(inner) = v.get("mpc_orb") {
        return parse_mpc_orb(inner);
    }
    let car = v.get("CAR").ok_or("no CAR (Cartesian) representation in the mpc_orb JSON")?;
    let vals: Vec<f64> = car["coefficient_values"].as_array().ok_or("CAR has no coefficient_values")?.iter().filter_map(|x| x.as_f64()).collect();
    if vals.len() < 6 {
        return Err("CAR has fewer than 6 values".into());
    }
    let names: Vec<String> = car["coefficient_names"].as_array().map(|a| a.iter().filter_map(|x| x.as_str().map(|s| s.to_string())).collect()).unwrap_or_default();
    let cov_obj = &car["covariance"];
    let mut cov_ecl = [0.0f64; 36];
    for i in 0..6 {
        for j in i..6 {
            let c = cov_obj.get(format!("cov{}{}", i, j)).and_then(|x| x.as_f64()).ok_or_else(|| format!("CAR covariance lacks cov{}{}", i, j))?;
            cov_ecl[i * 6 + j] = c;
            cov_ecl[j * 6 + i] = c;
        }
    }
    let mut nongrav = [0.0; 3];
    for (k, n) in ["A1", "A2", "A3"].iter().enumerate() {
        if let Some(pos) = names.iter().position(|s| s == n) {
            nongrav[k] = vals.get(pos).copied().unwrap_or(0.0);
        }
    }
    let ep = &v["epoch_data"];
    let epoch = ep["epoch"].as_f64().ok_or("no epoch in epoch_data")?;
    let jd = match ep["timeform"].as_str().unwrap_or("MJD").to_uppercase().as_str() {
        "MJD" => epoch + 2400000.5,
        _ => epoch,
    };
    let epoch_tdb = match ep["timesystem"].as_str().unwrap_or("TT").to_uppercase().as_str() {
        "TDB" => jd,
        "TT" | "TDT" => tt_to_tdb(jd),
        other => return Err(format!("unsupported time system {} in mpc_orb", other).into()),
    };
    let refsys = v["system_data"]["refsys"].as_str().unwrap_or("Ecliptic").to_lowercase();
    let m = if refsys.starts_with("eclip") { ecliptic_to_j2000() } else { Matrix3::identity() };
    let p = m * Vector3::new(vals[0], vals[1], vals[2]);
    let vv = m * Vector3::new(vals[3], vals[4], vals[5]);
    // Rotate the covariance: R C R^T, R = blockdiag(m, m).
    let mut r6 = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            r6[i][j] = m[(i, j)];
            r6[i + 3][j + 3] = m[(i, j)];
        }
    }
    let mut cov = [0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                for l in 0..6 {
                    s += r6[i][k] * cov_ecl[k * 6 + l] * r6[j][l];
                }
            }
            cov[i * 6 + j] = s;
        }
    }
    let dd = &v["designation_data"];
    let name = dd["permid"]
        .as_str()
        .or_else(|| dd["unpacked_primary_provisional_designation"].as_str())
        .or_else(|| dd["orbfit_name"].as_str())
        .unwrap_or("")
        .trim_matches(|c| c == '(' || c == ')')
        .to_string();
    let mag = &v["magnitude_data"];
    Ok(MpcOrb {
        name,
        epoch_tdb,
        helio_j2000: [p.x, p.y, p.z, vv.x, vv.y, vv.z],
        covariance: cov,
        h: mag["H"].as_f64().unwrap_or(f64::NAN),
        g: mag["G"].as_f64().unwrap_or(0.15),
        u: v["orbit_fit_statistics"]["U_param"].as_f64().unwrap_or(f64::NAN),
        nongrav,
    })
}

/// The in-orbit longitude runoff (radians per decade) for the MPC uncertainty parameter U: the
/// upper end of U's bin, `exp(1.49 U)` arcseconds (MPC: U = int(ln(runoff) / 1.49) + 1;
/// runoff < 1" is U = 0). U is derived from the uncertainties of the perihelion time and the
/// period, so for very long periods (Sedna: U = 5) it overstates the near-term uncertainty.
pub fn runoff_from_u(u: f64) -> f64 {
    let u = u.clamp(0.0, 9.0);
    (1.49 * u).exp() * std::f64::consts::PI / (180.0 * 3600.0)
}
