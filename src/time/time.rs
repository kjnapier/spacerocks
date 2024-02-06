use pyo3::prelude::*;
use pyo3::types::PyType;

use std::ops::{AddAssign, Add};
use std::collections::HashMap;
use chrono::{DateTime, TimeZone, Utc};
use crate::time::leapseconds::LEAP_SECONDS;
use lazy_static::lazy_static;


#[pyclass]
pub struct Time {
    #[pyo3(get, set)]
    pub epoch: f64,
    #[pyo3(get, set)]
    pub timescale: String,
    #[pyo3(get, set)]
    pub format: String,
}

#[pymethods]
impl Time {

    #[new]
    pub fn new(epoch: f64, timescale: &str, format: &str) -> Self {
        Time {
            epoch: epoch,
            timescale: timescale.to_string(),
            format: format.to_string(),
        }
    }
    
    #[classmethod]
    pub fn now(cls: &PyType) -> Self {
        let now = Utc::now();
        let x = now.to_rfc3339_opts(chrono::SecondsFormat::Millis, true);
        Time {
            epoch: isot_to_julian(&x),
            timescale: "utc".to_string(),
            format: "jd".to_string(),
        }
    }

    pub fn to_utc(&mut self) {
        if self.timescale == "utc" {
            return;
        }
        self.epoch = tdb_to_utc(self.epoch);
        self.timescale = "utc".to_string();
    }

    pub fn to_tdb(&mut self) {
        if self.timescale == "tdb" {
            return;
        }
        self.epoch = utc_to_tdb(self.epoch);
        self.timescale = "tdb".to_string();
    }

    pub fn change_timescale(&mut self, timescale: &str) {
        if timescale == "utc" {
            self.to_utc();
        } else if timescale == "tdb" {
            self.to_tdb();
        } else {
            panic!("Invalid timescale: {}", timescale);
        }
    }

    pub fn calendar(&self) -> String {
        return jd_to_calendar(&self.epoch);
    }
}

impl AddAssign<f64> for Time {
    fn add_assign(&mut self, dt: f64) {
        self.epoch += dt;
    }
}


fn utc_to_tdb(epoch: f64) -> f64 {
    // to tai
    let leapseconds = get_leap_seconds_at_epoch(epoch);
    let mut epoch = epoch;

    epoch += leapseconds / 86400.0;

    // to tt
    epoch += 32.184 / 86400.0;

    // to tdb
    let g = (357.53 + 0.9856003 * (epoch - 2451545.0)).to_radians();
    epoch += (0.001658 * g.sin() + 0.000014 * (2.0 * g).sin()) / 86400.0;
    return epoch;
}

fn tdb_to_utc(epoch: f64) -> f64 {
    // to  tt
    let g = (357.53 + 0.9856003 * (epoch - 2451545.0)).to_radians();
    let mut epoch = epoch;
    epoch -= (0.001658 * g.sin() + 0.000014 * (2.0 * g).sin()) / 86400.0;

    // to tai
    epoch -= 32.184 / 86400.0;

    // to utc
    let leapseconds = get_leap_seconds_at_epoch(epoch);
    epoch -= leapseconds / 86400.0;
    return epoch;
}


fn get_leap_seconds_at_epoch(jd: f64) -> f64 {
   
    let mut num_leap_seconds = 0.0;
    for &(time, leap_seconds) in &LEAP_SECONDS {
        if jd >= time {
            num_leap_seconds = leap_seconds;
            break;
        }
    }
    return num_leap_seconds;
}

fn get_isot_now() -> String {
    let now = Utc::now();
    let x = now.to_rfc3339_opts(chrono::SecondsFormat::Millis, true);
    return x;
}

fn isot_to_julian(isot: &str) -> f64 {
    let datetime: DateTime<Utc> = Utc.datetime_from_str(isot, "%Y-%m-%dT%H:%M:%S%.fZ").unwrap();
    let unix_time = datetime.timestamp() as f64;
    let julian_day = unix_time / 86400.0 + 2440587.5;
    julian_day
}


// hash mapping integers to month name
lazy_static! {
    static ref MONTHS: HashMap<u32, &'static str> = {
        let mut m = HashMap::new();
        m.insert(1, "Jan");
        m.insert(2, "Feb");
        m.insert(3, "Mar");
        m.insert(4, "Apr");
        m.insert(5, "May");
        m.insert(6, "Jun");
        m.insert(7, "Jul");
        m.insert(8, "Aug");
        m.insert(9, "Sep");
        m.insert(10, "Oct");
        m.insert(11, "Nov");
        m.insert(12, "Dec");
        m
    };
}

fn jd_to_calendar(jd: &f64) -> String {
    let jd = jd + 0.5;
    let z = jd.trunc() as i32;
    let a = if z < 2299161 {
        z
    } else {
        let alpha = ((z as f64 - 1867216.25) / 36524.25).floor() as i32;
        z + 1 + alpha - (alpha / 4)
    };
    let b = a + 1524;
    let c = ((b as f64 - 122.1) / 365.25).floor() as i32;
    let d = (365.25 * c as f64).floor() as i32;
    let e = ((b as f64 - d as f64) / 30.6001).floor() as u32;
    let day = b - d - ((30.6001 * e as f64) as i32);
    let month = if e < 14 {
        e - 1
    } else {
        e - 13
    };
    let year = if month > 2 {
        c - 4716
    } else {
        c - 4715
    };
    format!("{} {} {}", day, MONTHS.get(&month).unwrap(), year)
}