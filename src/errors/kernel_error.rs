#[derive(Debug, PartialEq)]
/// Errors that can occur when managing SPICE kernels.
pub enum KernelError {
    /// Failed to read or write file
    IoError(String),
    /// Failed to parse configuration
    ConfigError(String),
    /// Failed to download kernel
    DownloadError(String),
    /// Kernel not found in specified locations
    KernelNotFound(String),
    /// Invalid configuration provided
    InvalidConfig(String),
}

impl std::fmt::Display for KernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KernelError::IoError(msg) => write!(f, "IO error: {}", msg),
            KernelError::ConfigError(msg) => write!(f, "Configuration error: {}", msg),
            KernelError::DownloadError(msg) => write!(f, "Download error: {}", msg),
            KernelError::KernelNotFound(name) => write!(f, "Kernel not found: {}", name),
            KernelError::InvalidConfig(msg) => write!(f, "Invalid configuration: {}", msg),
        }
    }
}

impl std::error::Error for KernelError {}

