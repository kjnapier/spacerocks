use serde::{Serialize, Deserialize};


/// Represents different time formats for astronomical calculations
///
/// Supported formats:
/// - JD (Julian Date)
/// - MJD (Modified Julian Date)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[derive(Default)]
pub enum TimeFormat {
    #[default]
    JD,
    MJD,
}

impl TimeFormat {

    /// Returns a slice of all valid time format string representations
    ///
    /// Used for validating input strings and providing suggestions for invalid inputs
    pub fn variants() -> &'static [&'static str] {
        &["JD", "MJD"]
    }

    /// Creates a TimeFormat from a string representation
    ///
    /// # Arguments
    /// * `s` - String slice representing the time format
    ///
    /// # Returns
    /// * Some(TimeFormat) if the string is valid, None otherwise
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "JD" => Some(TimeFormat::JD),
            "MJD" => Some(TimeFormat::MJD),
            _ => None,
        }
    }

    /// Converts the time format to its string representation
    ///
    /// # Returns
    /// * A string slice representing the time format
    pub fn to_str(&self) -> &str {
        match self {
            TimeFormat::JD => "JD",
            TimeFormat::MJD => "MJD",
        }
    }
}

impl std::fmt::Display for TimeFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            TimeFormat::JD => write!(f, "JD"),
            TimeFormat::MJD => write!(f, "MJD"),
        }
    }
}

