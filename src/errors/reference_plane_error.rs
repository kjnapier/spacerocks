#[derive(Debug)]
/// Error from parsing an invalid reference plane string.
pub enum ReferencePlaneError {
    InvalidReferencePlane(String),
}

impl std::fmt::Display for ReferencePlaneError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReferencePlaneError::InvalidReferencePlane(s) => write!(f, "Invalid origin: {}", s),
        }
    }
}

impl std::error::Error for ReferencePlaneError {}