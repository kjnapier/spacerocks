//! Handling various time formats and scales.
//!
//! This module provides functionality for:
//! - Converting between different time scales (UTC, TDB, TT, TAI)
//! - Working with different time formats (JD, MJD)
//! - Handling leap seconds, UTC's rate offsets from 1960 to 1972, and ΔT (TT − UT) before 1960
//! - Converting between calendar dates and Julian dates
//!
//! # Example
//! ```rust
//! use spacerocks::time::Time;
//!
//! let time = Time::now();
//! println!("Current time: {}", time.calendar());
//! ```

pub mod time;
pub use self::time::{Time};

pub mod timeformat;
pub use self::timeformat::TimeFormat;

pub mod timescale;
pub use self::timescale::TimeScale;

pub mod leapseconds;

pub mod deltat;
pub use self::deltat::delta_t;

pub mod conversions;
pub use self::conversions::*; 