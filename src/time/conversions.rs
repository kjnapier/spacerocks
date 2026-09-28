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

/// 1972 January 1, 0h UTC (JD): from here on, TAI − UTC is a whole number of leap seconds.
pub const UTC_LEAP_SECOND_ERA: f64 = 2441317.5;

/// 1960 January 1, 0h UTC (JD): the start of UTC, as SOFA/ERFA's `dat` has it. Before this
/// date "UTC" is taken to be UT, and TT − UT = ΔT (see [`crate::time::deltat`]).
pub const UTC_START: f64 = 2436934.5;

/// TAI − UTC from 1960 to 1972, when UTC ran at an offset rate from TAI (SOFA/ERFA `dat`, from
/// the USNO table): from each JD on, TAI − UTC = offset + (MJD − reference MJD) × rate seconds.
const UTC_1960_1972: [(f64, f64, f64, f64); 14] = [
    (2436934.5, 1.4178180, 37300.0, 0.0012960), // 1960 Jan 1
    (2437300.5, 1.4228180, 37300.0, 0.0012960), // 1961 Jan 1
    (2437512.5, 1.3728180, 37300.0, 0.0012960), // 1961 Aug 1
    (2437665.5, 1.8458580, 37665.0, 0.0011232), // 1962 Jan 1
    (2438334.5, 1.9458580, 37665.0, 0.0011232), // 1963 Nov 1
    (2438395.5, 3.2401300, 38761.0, 0.0012960), // 1964 Jan 1
    (2438486.5, 3.3401300, 38761.0, 0.0012960), // 1964 Apr 1
    (2438639.5, 3.4401300, 38761.0, 0.0012960), // 1964 Sep 1
    (2438761.5, 3.5401300, 38761.0, 0.0012960), // 1965 Jan 1
    (2438820.5, 3.6401300, 38761.0, 0.0012960), // 1965 Mar 1
    (2438942.5, 3.7401300, 38761.0, 0.0012960), // 1965 Jul 1
    (2439004.5, 3.8401300, 38761.0, 0.0012960), // 1965 Sep 1
    (2439126.5, 4.3131700, 39126.0, 0.0025920), // 1966 Jan 1
    (2439887.5, 4.2131700, 39126.0, 0.0025920), // 1968 Feb 1
];

/// TT − UTC in seconds at the UTC Julian date `jd`, before 1972: from 1960, 32.184 s plus
/// TAI − UTC from the 1960–1972 table; before 1960, when there was no UTC and times are UT,
/// ΔT (Stephenson, Morrison & Hohenkerk).
fn tt_minus_utc_before_1972(jd: f64) -> f64 {
    if jd < UTC_START {
        return crate::time::deltat::delta_t(jd);
    }
    let i = UTC_1960_1972.partition_point(|r| r.0 <= jd) - 1;
    let (_, offset, mjd_ref, rate) = UTC_1960_1972[i];
    32.184 + offset + (jd - 2400000.5 - mjd_ref) * rate
}

/// TT − UTC in seconds at the UTC Julian date `jd`: 32.184 s plus the leap seconds from 1972 on,
/// the 1960–1972 rate offsets before that, and ΔT = TT − UT before 1960.
pub fn tt_minus_utc(jd: f64) -> f64 {
    if jd >= UTC_LEAP_SECOND_ERA {
        32.184 + get_leap_seconds_at_epoch(jd)
    } else {
        tt_minus_utc_before_1972(jd)
    }
}

/// Converts UTC (Universal Time Coordinated) to TAI (International Atomic Time)
///
/// From 1972, TAI − UTC is the leap-second count. From 1960 to 1972 it follows UTC's rate
/// offsets (as SOFA/ERFA). Before 1960, when there was no UTC, the epoch is taken as UT, and
/// TAI = TT − 32.184 s with TT − UT = ΔT.
///
/// # Arguments
/// 
/// * `epoch` - The epoch in UTC Julian Date
/// 
/// # Returns
/// 
/// * The epoch in TAI Julian Date
pub fn utc_to_tai(epoch: f64) -> f64 {
    if epoch >= UTC_LEAP_SECOND_ERA {
        let leapseconds = get_leap_seconds_at_epoch(epoch);
        epoch + leapseconds / 86400.0
    } else {
        epoch + (tt_minus_utc_before_1972(epoch) - 32.184) / 86400.0
    }
}

/// Converts TAI (International Atomic Time) to UTC (Universal Time Coordinated): the inverse
/// of [`utc_to_tai`].
/// 
/// # Arguments
/// 
/// * `epoch` - The epoch in TAI Julian Date
/// 
/// # Returns
/// 
/// * The epoch in UTC Julian Date
pub fn tai_to_utc(epoch: f64) -> f64 {
    if epoch >= UTC_LEAP_SECOND_ERA + 10.0 / 86400.0 {
        let leapseconds = get_leap_seconds_at_epoch(epoch);
        return epoch - leapseconds / 86400.0;
    }
    // TAI − UTC changes by at most ~1e-7 s per second here, so this converges at once.
    let mut utc = epoch;
    for _ in 0..4 {
        utc = epoch - (tt_minus_utc_before_1972(utc) - 32.184) / 86400.0;
    }
    utc
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
/// # use spacerocks::time::*;
/// let jd = 2451545.0;
/// let date = jd_to_calendar(&jd);
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
/// # use spacerocks::time::*;
/// let julian_date = isot_to_julian("2024-12-11T12:34:56.789Z");
/// println!("Julian Date: {}", julian_date);
/// ```
pub fn isot_to_julian(isot: &str) -> f64 {
    let datetime: DateTime<Utc> = chrono::NaiveDateTime::parse_from_str(isot.trim_end_matches('Z'), "%Y-%m-%dT%H:%M:%S%.f")
        .unwrap()
        .and_utc();
    // Keep the fractional seconds (timestamp() alone truncates them).
    let unix_time = datetime.timestamp() as f64 + datetime.timestamp_subsec_nanos() as f64 * 1e-9;
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tt_minus_utc_by_era() {
        // leap seconds (unchanged)
        assert_eq!(tt_minus_utc(2451545.0), 32.184 + 32.0);
        assert_eq!(tt_minus_utc(UTC_LEAP_SECOND_ERA), 42.184);
        // 1960-1972 rate offsets (ERFA dat): 1965 Jan 1, and the day before the 1972 step
        assert!((tt_minus_utc(2438761.5) - (32.184 + 3.5401300)).abs() < 1e-12);
        let dec31 = 2441316.5;
        assert!((tt_minus_utc(dec31) - (32.184 + 4.2131700 + (41316.0 - 39126.0) * 0.002592)).abs() < 1e-9);
        // steps apply at 0h: 1961 Aug 1 lowered TAI - UTC by 0.05 s
        assert!((tt_minus_utc(2437512.5 - 1e-6) - tt_minus_utc(2437512.5) - 0.05).abs() < 1e-6);
        // before 1960: Delta T
        assert_eq!(tt_minus_utc(2415020.5), crate::time::deltat::delta_t(2415020.5));
    }

    #[test]
    fn utc_tai_round_trips_in_every_era() {
        for jd in [2305447.5, 2415020.3, 2436934.2, 2437000.7, 2438500.1, 2441000.9, 2451545.0, 2460000.25] {
            let back = tai_to_utc(utc_to_tai(jd));
            assert!((back - jd).abs() * 86400.0 < 1e-4, "{}: {} s", jd, (back - jd) * 86400.0);
            let expected = tt_minus_utc(jd) - 32.184;
            assert!(((utc_to_tai(jd) - jd) * 86400.0 - expected).abs() < 1e-4);
        }
    }
}
