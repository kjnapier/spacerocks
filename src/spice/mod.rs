//! SPICE kernel management module.
//! 
//! Provides functionality for managing SPICE kernels including:
//! - Configuration handling via TOML files
//! - Kernel loading and unloading
//! - Automatic downloading of missing kernels

// pub mod spicekernel;
// pub use self::spicekernel::SpiceKernel;

pub mod config;
pub mod spicekernel;

pub use spicekernel::SpiceKernel;
pub use config::KernelSpec;