#[derive(Debug)]
/// Errors that can occur when working with time objects.
pub enum TimeError {
    InvalidTimeScale(String),
    InvalidTimeFormat(String),
}

impl std::fmt::Display for TimeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimeError::InvalidTimeScale(s) => write!(f, "Invalid timescale: {}. Needs to be 'utc', 'tdb', 'tt', or 'tai'.", s),
            TimeError::InvalidTimeFormat(s) => write!(f, "Invalid time format: {}. Needs to be 'jd' or 'mjd'.", s),
        }
    }
}

impl std::error::Error for TimeError {}