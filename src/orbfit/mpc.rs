//! Observations from the Minor Planet Center: fetch an object's optical astrometry (or a
//! random numbered asteroid's) and turn it into [`Astrometry`].
//!
//! [`fetch_observations`] asks the MPC observations API for an object's 80-column (Obs80)
//! records, [`parse_obs80`] reads them, and [`mpc_astrometry`] keeps the ground-based optical
//! ones in J2000, converts their times to TDB, weights them by [`veres_sigma`] and places the
//! stations with the kernel. [`fetch_astrometry`] does all three, for a designation or for a
//! random numbered asteroid ([`random_numbered`]).
use rand::Rng;
use rayon::prelude::*;
use serde_json::{json, Value};

use crate::orbfit::{observer_positions, veres_sigma, Astrometry, ARCSEC};
use crate::spice::SpiceKernel;
use crate::time::Time;
use crate::OBSERVATORIES;

/// The MPC observations API (`GET` with a JSON body).
pub const MPC_OBSERVATIONS_API: &str = "https://data.minorplanetcenter.net/api/get-obs";

/// Numbered asteroids [`random_numbered`] picks from (all of them are numbered well past this).
pub const RANDOM_NUMBERED_MAX: u32 = 750_000;

/// Fewest usable detections a random object needs before [`fetch_astrometry`] takes it.
const RANDOM_MIN_DETECTIONS: usize = 20;

/// Random picks [`fetch_astrometry`] tries before giving up.
const RANDOM_TRIES: usize = 5;

/// Obs80 note-2 codes [`parse_obs80`] skips: satellite, roving and radar records (which take a
/// second line or aren't positions), deleted and replaced discoveries, occultations, offsets,
/// and B1950 positions (`A`).
const SKIPPED_NOTES: &str = "SsVvRrXxEOA";

/// Station codes that need more than the Obs80 line (roving, space-based, geocentric).
const SPECIAL_STATIONS: [&str; 9] = ["244", "245", "247", "248", "249", "250", "270", "275", "500"];

/// One optical observation from an Obs80 record.
#[derive(Debug, Clone, PartialEq)]
pub struct MpcObservation {
    /// UTC Julian date.
    pub jd_utc: f64,
    /// J2000 right ascension and declination (radians).
    pub ra: f64,
    pub dec: f64,
    /// MPC station code.
    pub station: String,
    /// Star catalog code (column 72), empty if none.
    pub catalog: String,
    /// Program code (column 14), empty if none.
    pub program: String,
}

/// Julian date of a Gregorian calendar date (day with its fraction).
fn julian_date(year: i32, month: i32, day: f64) -> f64 {
    let (y, m) = if month <= 2 { (year - 1, month + 12) } else { (year, month) };
    let a = y.div_euclid(100);
    let b = 2 - a + a.div_euclid(4);
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (m as f64 + 1.0)).floor() + day + b as f64 - 1524.5
}

/// Three space-separated sexagesimal fields, like `12 34 56.789`, as `(sign, a, b, c)`.
fn sexagesimal(s: &str) -> Option<(f64, f64, f64, f64)> {
    let s = s.trim();
    let (sign, s) = match s.as_bytes().first()? {
        b'-' => (-1.0, &s[1..]),
        b'+' => (1.0, &s[1..]),
        _ => (1.0, s),
    };
    let mut f = s.split_whitespace().map(|v| v.parse::<f64>().ok());
    let (a, b, c) = (f.next()??, f.next()??, f.next()??);
    Some((sign, a, b, c))
}

/// Parse one Obs80 optical record, or `None` for a line this module doesn't use.
fn parse_obs80_line(line: &str) -> Option<MpcObservation> {
    let line = line.trim_end_matches(['\r', '\n']);
    if line.len() < 80 || !line.is_ascii() {
        return None;
    }
    let col = |a: usize, b: usize| &line[a - 1..b];
    let note2 = col(15, 15).chars().next()?;
    if SKIPPED_NOTES.contains(note2) {
        return None;
    }
    let station = col(78, 80).trim().to_string();
    if SPECIAL_STATIONS.contains(&station.as_str()) || OBSERVATORIES.get(station.as_str()).is_none() {
        return None;
    }
    let mut date = col(16, 32).split_whitespace();
    let year: i32 = date.next()?.parse().ok()?;
    let month: i32 = date.next()?.parse().ok()?;
    let day: f64 = date.next()?.parse().ok()?;
    let (_, h, m, s) = sexagesimal(col(33, 44))?;
    let (sign, d, dm, ds) = sexagesimal(col(45, 56))?;
    let ra = (15.0 * (h + m / 60.0 + s / 3600.0)).to_radians();
    let dec = (sign * (d + dm / 60.0 + ds / 3600.0)).to_radians();
    if !(0.0..2.0 * std::f64::consts::PI).contains(&ra) || dec.abs() > std::f64::consts::FRAC_PI_2 || !(1..=12).contains(&month) {
        return None;
    }
    Some(MpcObservation {
        jd_utc: julian_date(year, month, day),
        ra,
        dec,
        station,
        catalog: col(72, 72).trim().to_string(),
        program: col(14, 14).trim().to_string(),
    })
}

/// The ground-based optical J2000 observations among Obs80 records (one per line). Satellite,
/// roving, radar, occultation, offset, deleted and B1950 records, and stations not in the
/// observatory table, are left out.
pub fn parse_obs80(text: &str) -> Vec<MpcObservation> {
    text.lines().filter_map(parse_obs80_line).collect()
}

/// The Obs80 text in an MPC observations API response: `[{"OBS80": "..."}]`, where the text may
/// also come as `{"data": "..."}` or a list of lines.
pub fn obs80_from_response(response: &Value) -> Result<String, String> {
    let item = match response {
        Value::Array(items) => items.first().ok_or("the MPC has no observations of that object")?,
        other => other,
    };
    if let Some(e) = item.get("error").or_else(|| item.get("message")).and_then(Value::as_str) {
        return Err(format!("MPC: {e}"));
    }
    let obs = item.get("OBS80").or_else(|| item.get("obs80")).ok_or("the MPC response has no OBS80 field")?;
    let text = match obs {
        Value::String(s) => s.clone(),
        Value::Object(o) => o.get("data").and_then(Value::as_str).ok_or("the MPC response's OBS80 has no data")?.to_string(),
        Value::Array(lines) => lines.iter().filter_map(Value::as_str).collect::<Vec<_>>().join("\n"),
        _ => return Err("the MPC response's OBS80 isn't text".into()),
    };
    Ok(text)
}

/// Fetch the Obs80 records of the object `designation` (a number, a name or a provisional
/// designation) from the MPC observations API.
pub fn fetch_obs80(designation: &str) -> Result<String, String> {
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(60))
        .user_agent(concat!("spacerocks/", env!("CARGO_PKG_VERSION")))
        .build()
        .map_err(|e| e.to_string())?;
    let body = json!({ "desigs": [designation], "output_format": ["OBS80"] });
    let response = client.get(MPC_OBSERVATIONS_API).json(&body).send().map_err(|e| format!("couldn't reach the MPC: {e}"))?;
    let status = response.status();
    let text = response.text().map_err(|e| format!("couldn't read the MPC's answer: {e}"))?;
    if !status.is_success() {
        let detail: String = text.chars().take(200).collect();
        return Err(format!("MPC answered {status} for '{designation}': {}", detail.trim()));
    }
    let value: Value = serde_json::from_str(&text).map_err(|e| format!("the MPC's answer isn't JSON: {e}"))?;
    obs80_from_response(&value)
}

/// The ground-based optical observations of `designation`, from the MPC.
pub fn fetch_observations(designation: &str) -> Result<Vec<MpcObservation>, String> {
    Ok(parse_obs80(&fetch_obs80(designation)?))
}

/// A random numbered asteroid, as its number.
pub fn random_numbered() -> String {
    rand::thread_rng().gen_range(1..=RANDOM_NUMBERED_MAX).to_string()
}

/// Astrometry from MPC observations: times to TDB, Vereš et al. (2017) uncertainties, station
/// positions from `kernel`. Observations the kernel can't place (outside its time span) are
/// dropped. With more than `max_detections`, keeps that many spread evenly over the arc
/// (the first and last included).
pub fn mpc_astrometry(observations: &[MpcObservation], kernel: &SpiceKernel, max_detections: usize) -> Result<Astrometry, String> {
    let mut obs: Vec<&MpcObservation> = observations.iter().collect();
    obs.sort_by(|a, b| a.jd_utc.total_cmp(&b.jd_utc));
    obs.dedup_by(|a, b| a.jd_utc == b.jd_utc && a.station == b.station);
    let placed: Vec<Result<(&MpcObservation, f64, [f64; 3]), String>> = obs
        .par_iter()
        .map(|o| {
            let jd = Time::new(o.jd_utc, "utc", "jd").map_err(|e| e.to_string())?.tdb().jd();
            let p = observer_positions(&[o.station.clone()], &[jd], kernel)?;
            Ok((*o, jd, p[0]))
        })
        .collect();
    let first_error = placed.iter().find_map(|p| p.as_ref().err().cloned());
    let placed: Vec<_> = placed.into_iter().filter_map(Result::ok).filter(|(_, _, p)| p.iter().all(|v| v.is_finite())).collect();
    if placed.is_empty() {
        return Err(match first_error {
            Some(e) => format!("the kernel can't place the observing stations: {e}"),
            None => "no usable ground-based optical observations".into(),
        });
    }
    let keep: Vec<usize> = if max_detections >= 2 && placed.len() > max_detections {
        let step = (placed.len() - 1) as f64 / (max_detections - 1) as f64;
        let mut k: Vec<usize> = (0..max_detections).map(|i| (i as f64 * step).round() as usize).collect();
        k.dedup();
        k
    } else {
        (0..placed.len()).collect()
    };
    let (mut epoch, mut ra, mut dec, mut sigma, mut observer) = (vec![], vec![], vec![], vec![], vec![]);
    for &i in &keep {
        let (o, jd, p) = placed[i];
        epoch.push(jd);
        ra.push(o.ra);
        dec.push(o.dec);
        sigma.push(veres_sigma(&o.station, jd, Some(&o.catalog), Some(&o.program)) * ARCSEC);
        observer.push(p);
    }
    Astrometry::new(epoch, ra, dec, observer)?.with_sigma(&sigma, &sigma)
}

/// Fetch the astrometry of `designation` from the MPC, or with `None` of a random numbered
/// asteroid (trying another if one has fewer than 20 usable detections). Returns the object's
/// label and its astrometry (see [`mpc_astrometry`]).
pub fn fetch_astrometry(designation: Option<&str>, kernel: &SpiceKernel, max_detections: usize) -> Result<(String, Astrometry), String> {
    match designation.map(str::trim).filter(|d| !d.is_empty()) {
        Some(d) => Ok((d.to_string(), mpc_astrometry(&fetch_observations(d)?, kernel, max_detections)?)),
        None => {
            let mut last = String::new();
            for _ in 0..RANDOM_TRIES {
                let number = random_numbered();
                match fetch_observations(&number).and_then(|o| mpc_astrometry(&o, kernel, max_detections)) {
                    Ok(a) if a.len() >= RANDOM_MIN_DETECTIONS => return Ok((format!("({number})"), a)),
                    Ok(a) => last = format!("({number}) has only {} usable detections", a.len()),
                    // Unreachable or refused: another number won't help.
                    Err(e) if e.starts_with("couldn't reach") => return Err(e),
                    Err(e) => last = e,
                }
            }
            Err(format!("no usable random object in {RANDOM_TRIES} tries ({last})"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn julian_dates() {
        assert_eq!(julian_date(2000, 1, 1.5), 2451545.0);
        assert_eq!(julian_date(1987, 6, 19.5), 2446966.0);
        assert_eq!(julian_date(1600, 12, 31.0), 2305812.5);
    }
}
