//! Error types for the spacerocks library.
//!
//! This module provides a comprehensive set of error types for
//! handling various failure conditions that may arise when using Spacerocks.
//!
//! Each error type is specialized for a particular domain:
//! - [`OrbitError`]: Orbital mechanics calculations
//! - [`OriginError`]: Coordinate system origins
//! - [`TimeError`]: Time-scale and epoch handling
//! - [`SimulationError`]: N-body simulation operations
//! - [`ReferencePlaneError`]: Reference frame transformations
//! - [`KernelError`]: SPICE kernel management
//!
//! # Examples
//!
//! ```
//! use spacerocks::{OrbitError, TimeError};
//!
//! fn propagate_orbit(epoch: f64) -> Result<(), Box<dyn std::error::Error>> {
//!     // Example error handling for orbital calculations
//!     if epoch < 0.0 {
//!         return Err(Box::new(TimeError::InvalidEpoch(epoch)));
//!     }
//!     Ok(())
//! }
//! ```

pub mod orbit_error;
pub use self::orbit_error::OrbitError;

pub mod origin_error;
pub use self::origin_error::OriginError;

pub mod time_error;
pub use self::time_error::TimeError;

pub mod simulation_error;
pub use self::simulation_error::SimulationError;

pub mod reference_plane_error;
pub use self::reference_plane_error::ReferencePlaneError;

pub mod kernel_error;
pub use self::kernel_error::KernelError;
