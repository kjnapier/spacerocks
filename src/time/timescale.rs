use serde::{Serialize, Deserialize};

/// Represents different time scales used in astronomical calculations
///
/// Supported time scales:
/// - UTC (Universal Time Coordinated)
/// - TDB (Barycentric Dynamical Time)
/// - TT (Terrestrial Time)
/// - TAI (International Atomic Time)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[derive(Default)]
pub enum TimeScale {
    #[default]
    UTC,
    TDB,
    TT, 
    TAI,
}

impl TimeScale {

    /// Returns a slice of all valid time scale string representations
    ///
    /// Used for validating input strings and providing suggestions for invalid inputs
    pub fn variants() -> &'static [&'static str] {
        &["UTC", "TDB", "TT", "TAI"] 
    }

    /// Converts the time scale to its string representation
    ///
    /// # Returns
    /// * A string slice representing the time scale
    pub fn to_str(&self) -> &str {
        match self {
            TimeScale::UTC => "UTC",
            TimeScale::TDB => "TDB",
            TimeScale::TT => "TT",
            TimeScale::TAI => "TAI",
        }
    }
}

impl std::fmt::Display for TimeScale {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            TimeScale::UTC => write!(f, "UTC"),
            TimeScale::TDB => write!(f, "TDB"),
            TimeScale::TT => write!(f, "TT"),
            TimeScale::TAI => write!(f, "TAI"),
        }
    }
}

