#[derive(Debug)]
/// Error from parsing an invalid origin string.
pub enum OriginError {
    InvalidOrigin(String),
}

impl std::fmt::Display for OriginError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OriginError::InvalidOrigin(s) => write!(f, "Invalid origin: {}", s),
        }
    }
}

impl std::error::Error for OriginError {}