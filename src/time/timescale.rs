use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum TimeScale {
    UTC,
    TDB,
}

impl TimeScale {
    pub fn to_str(&self) -> &str {
        match self {
            TimeScale::UTC => "UTC",
            TimeScale::TDB => "TDB",
        }
    }
}

impl std::fmt::Display for TimeScale {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            TimeScale::UTC => write!(f, "UTC"),
            TimeScale::TDB => write!(f, "TDB"),
        }
    }
}

impl Default for TimeScale {
    fn default() -> Self {
        TimeScale::UTC
    }
}
