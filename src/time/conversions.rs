//! Time scale conversion functions
//!
//! This module provides functions for converting between different astronomical time scales:
//! - UTC (Universal Time Coordinated)
//! - TAI (International Atomic Time)
//! - TT (Terrestrial Time)
//! - TDB (Barycentric Dynamical Time)
//!
//! All functions operate on Julian Dates.

use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::time::leapseconds::LEAP_SECONDS;
use chrono::{Utc, DateTime};
use chrono::TimeZone;

/// Converts UTC (Universal Time Coordinated) to TAI (International Atomic Time)
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in UTC Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TAI Julian Date
pub fn utc_to_tai(epoch: f64) -> f64 {
    let leapseconds = get_leap_seconds_at_epoch(epoch);
    epoch + leapseconds / 86400.0
}

/// Converts TAI (International Atomic Time) to UTC (Universal Time Coordinated)
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TAI Julian Date
/// 
/// # Returns
/// 
/// * The epoch in UTC Julian Date
pub fn tai_to_utc(epoch: f64) -> f64 {
    let leapseconds = get_leap_seconds_at_epoch(epoch);
    epoch - leapseconds / 86400.0
}

/// Converts TAI (International Atomic Time) to TT (Terrestrial Time)
/// TT differs from TAI by a constant offset of 32.184 seconds
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TAI Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TT Julian Date
pub fn tai_to_tt(epoch: f64) -> f64 {
    epoch + 32.184 / 86400.0
}

/// Converts TT (Terrestrial Time) to TAI (International Atomic Time)
/// TAI differs from TT by a constant offset of -32.184 seconds
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TT Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TAI Julian Date
pub fn tt_to_tai(epoch: f64) -> f64 {
    epoch - 32.184 / 86400.0
}

/// Converts TT (Terrestrial Time) to TDB (Barycentric Dynamical Time)
/// Includes periodic relativistic corrections
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TT Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TDB Julian Date
pub fn tt_to_tdb(epoch: f64) -> f64 {
    let g = (357.53 + 0.9856003 * (epoch - 2451545.0)).to_radians();
    epoch + (0.001658 * g.sin() + 0.000014 * (2.0 * g).sin()) / 86400.0
}

/// Converts TDB (Barycentric Dynamical Time) to TT (Terrestrial Time)
/// Removes periodic relativistic corrections
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TDB Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TT Julian Date
pub fn tdb_to_tt(epoch: f64) -> f64 {
    let g = (357.53 + 0.9856003 * (epoch - 2451545.0)).to_radians();
    epoch - (0.001658 * g.sin() + 0.000014 * (2.0 * g).sin()) / 86400.0
}

/// Converts UTC (Universal Time Coordinated) to TDB (Barycentric Dynamical Time)
/// Conversion chain: UTC -> TAI -> TT -> TDB
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in UTC Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TDB Julian Date
pub fn utc_to_tdb(epoch: f64) -> f64 {
    let tai = utc_to_tai(epoch);
    let tt = tai_to_tt(tai);
    tt_to_tdb(tt)
}

/// Converts TDB (Barycentric Dynamical Time) to UTC (Universal Time Coordinated)
/// Conversion chain: TDB -> TT -> TAI -> UTC
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TDB Julian Date
/// 
/// # Returns
/// 
/// * The epoch in UTC Julian Date
pub fn tdb_to_utc(epoch: f64) -> f64 {
    let tt = tdb_to_tt(epoch);
    let tai = tt_to_tai(tt);
    tai_to_utc(tai)
}

/// Converts UTC (Universal Time Coordinated) to TT (Terrestrial Time)
/// Conversion chain: UTC -> TAI -> TT
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in UTC Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TT Julian Date
pub fn utc_to_tt(epoch: f64) -> f64 {
    let tai = utc_to_tai(epoch);
    tai_to_tt(tai)
}

/// Converts TT (Terrestrial Time) to UTC (Universal Time Coordinated)
/// Conversion chain: TT -> TAI -> UTC
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TT Julian Date
/// 
/// # Returns
/// 
/// * The epoch in UTC Julian Date
pub fn tt_to_utc(epoch: f64) -> f64 {
    let tai = tt_to_tai(epoch);
    tai_to_utc(tai)
}

// Calendar related conversions
//
// The following functions convert between Julian Date (JD) and the Gregorian calendar.

// /hash mapping integers to month name
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

/// Converts a Julian Date to the Gregorian calendar
/// 
/// # Arguments
///
/// * `jd` - The Julian Date
///
/// # Returns
///
/// * A string representing the date in the format "DD Mon YYYY"
/// 
/// # Example
///
/// ```
/// let jd = 2451545.0;
/// let date = jd_to_calendar(jd);
/// println!("Date: {}", date);
/// ```
pub fn jd_to_calendar(jd: &f64) -> String {
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

/// Get number of leap seconds at a given epoch
///
/// # Arguments
///
/// * `jd` - The Julian Date
///
/// # Returns
///
/// * The number of leap seconds at the given epoch
pub fn get_leap_seconds_at_epoch(jd: f64) -> f64 {
   
    let mut num_leap_seconds = 0.0;
    for &(time, leap_seconds) in &LEAP_SECONDS {
        if jd >= time {
            num_leap_seconds = leap_seconds;
            break;
        }
    }
    num_leap_seconds
}

/// Converts an ISO 8601 formatted timestamp to a Julian Date.
///
/// # Arguments
///
/// * `isot` - A string representing the ISO 8601 formatted timestamp 
///            (e.g., "2024-12-11T12:34:56.789Z") in UTC.
///
/// # Returns
///
/// * The Julian Date corresponding to the given timestamp.
///
/// # Example
///
/// ```
/// let julian_date = isot_to_julian("2024-12-11T12:34:56.789Z");
/// println!("Julian Date: {}", julian_date);
/// ```
pub fn isot_to_julian(isot: &str) -> f64 {
    let datetime: DateTime<Utc> = Utc.datetime_from_str(isot, "%Y-%m-%dT%H:%M:%S%.fZ").unwrap();
    // let datetime: DateTime<Utc> = DateTime::parse_from_str(isot, "%Y-%m-%dT%H:%M:%S%.fZ").unwrap().into();
    // let datetime: DateTime<Utc> = DateTime::parse_from_str(isot, "%Y-%m-%dT%H:%M:%S%.fZ").unwrap().into();
    let unix_time = datetime.timestamp() as f64;
    let julian_day = unix_time / 86400.0 + 2440587.5;
    julian_day
}




// use std::time::{SystemTime, UNIX_EPOCH};

// fn main() {
//     // Get the system's current time
//     let now = SystemTime::now();

//     // Calculate the duration since 1970-01-01 00:00:00 UTC
//     let duration_since_epoch = now
//         .duration_since(UNIX_EPOCH)
//         .expect("Time went backwards?");

//     // Convert that duration to whole seconds
//     let unix_timestamp = duration_since_epoch.as_secs();

//     println!("Current Unix timestamp: {}", unix_timestamp);
// }

