//! SPICE kernel management module.
//! 
//! Provides functionality for managing SPICE kernels including:
//! - Configuration handling via TOML files
//! - Kernel loading and unloading
//! - Automatic downloading of missing kernels

// pub mod spicekernel;
// pub use self::spicekernel::SpiceKernel;

pub mod config;

pub use spicekernel::SpiceKernel;
pub use config::KernelSpec;

pub mod spicebody;
pub mod spicekernel;
pub mod spk;
pub mod error;
pub mod bpc;

pub use spicebody::SpiceBody;
pub use spk::Spk;
pub use bpc::Bpc;
pub use error::SpiceError;
