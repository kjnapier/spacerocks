
use std::io;
use std::fmt;

use thiserror::Error;

/// Possible errors during SPK processing.
#[derive(Debug, Error)]
pub enum SpiceError {

    #[error("AstFile error")]
    AstFile,

    #[error("Body: {0} not found")]
    BodyNotFound(String),

    #[error("Nast error")]
    Nast,

    #[error("Coverage error")]
    Coverage,

    #[error("IO error: {0}")]
    IoError(io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),
}

impl From<io::Error> for SpiceError {
    fn from(err: io::Error) -> Self {
        SpiceError::IoError(err)
    }
}

// impl std::error::Error for SpiceError {}

// // Implement the Display trait for SpiceError.
// impl fmt::Display for SpiceError {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         match self {
//             SpiceError::AstFile => write!(f, "AstFile error"),
//             SpiceError::Nast => write!(f, "Nast error"),
//             SpiceError::Coverage => write!(f, "Coverage error"),
//             SpiceError::IoError(err) => write!(f, "IO error: {}", err),
//             SpiceError::ParseError(err) => write!(f, "Parse error: {}", err),
//         }
//     }
// }