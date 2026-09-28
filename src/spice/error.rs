//! Error type for the SPICE subsystem.

use std::io;
use thiserror::Error;

/// Errors produced while loading kernels or evaluating ephemerides/orientations.
#[derive(Debug, Error)]
pub enum SpiceError {
    /// Underlying I/O failure (file missing, unreadable, ...).
    #[error("I/O error for '{path}': {source}")]
    Io {
        path: String,
        #[source]
        source: io::Error,
    },

    /// The file is not a kernel we understand, or it is malformed.
    #[error("invalid kernel file '{path}': {reason}")]
    InvalidFile { path: String, reason: String },

    /// A segment type that is not implemented.
    #[error("{kind} segment type {data_type} is not supported (target/body {body})")]
    UnsupportedSegmentType {
        kind: &'static str,
        data_type: i32,
        body: i32,
    },

    /// A segment's contents are inconsistent.
    #[error("malformed {kind} type {data_type} segment for body {body}: {reason}")]
    MalformedSegment {
        kind: &'static str,
        data_type: i32,
        body: i32,
        reason: String,
    },

    /// Could not connect target and observer through the loaded SPK data at this epoch.
    #[error(
        "insufficient ephemeris data to compute the state of {target} relative to {observer} at ET {et} \
         (TDB seconds past J2000){detail}"
    )]
    InsufficientEphemerisData {
        target: i32,
        observer: i32,
        et: f64,
        detail: String,
    },

    /// No orientation data for a frame at this epoch.
    #[error("insufficient orientation data for frame '{frame}' at ET {et}: {detail}")]
    InsufficientOrientationData {
        frame: String,
        et: f64,
        detail: String,
    },

    /// A body name that could not be mapped to a NAIF ID.
    #[error("unknown body '{0}'")]
    UnknownBody(String),

    /// A frame name or ID that is not known.
    #[error("unknown frame '{0}'")]
    UnknownFrame(String),

    /// A frame whose class we do not implement (CK, dynamic, switch).
    #[error("frame '{frame}' has class {class}, which is not supported")]
    UnsupportedFrameClass { frame: String, class: i32 },

    /// Frame chain exceeded the maximum depth (probably a cycle in frame definitions).
    #[error("frame chain for '{0}' is too deep (cyclic frame definitions?)")]
    FrameChainTooDeep(String),

    /// Text kernel parsing problem.
    #[error("text kernel parse error in '{path}' line {line}: {reason}")]
    TextKernel {
        path: String,
        line: usize,
        reason: String,
    },

    /// A kernel in a kernel set could not be found locally (and was not downloaded).
    #[error("kernel '{name}' not found: {detail}")]
    KernelNotFound { name: String, detail: String },

    /// Downloading a kernel failed.
    #[error("failed to download {url}: {reason}")]
    Download { url: String, reason: String },

    /// Invalid kernel-set configuration.
    #[error("invalid kernel configuration: {0}")]
    Config(String),

    /// A required kernel pool variable is missing or has the wrong type.
    #[error("kernel pool variable '{name}': {reason}")]
    PoolVariable { name: String, reason: String },
}

impl SpiceError {
    pub(crate) fn io(path: &str, source: io::Error) -> Self {
        SpiceError::Io {
            path: path.to_string(),
            source,
        }
    }

    pub(crate) fn invalid(path: &str, reason: impl Into<String>) -> Self {
        SpiceError::InvalidFile {
            path: path.to_string(),
            reason: reason.into(),
        }
    }
}

pub type Result<T> = std::result::Result<T, SpiceError>;
