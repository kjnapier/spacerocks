//! Constants and observatory data.
//!
//! This module provides:
//! - Physical and astronomical constants in the `constants` module
//! - Observatory locations in the `observatories` module
//!
//! The constants include:
//! - Conversion factors (KM_TO_AU, M_TO_AU)
//! - Physical constants (SPEED_OF_LIGHT, GRAVITATIONAL_CONSTANT)
//! - Reference frame rotation matrices
//! - Solar system body masses
//!
//! Observatory locations are stored as three-component vectors where:
//! - First component is the observatory's longitude in radians
//! - Second and third components are direction cosines representing the
//!   observatory's position relative to Earth's center
//!
//! # Examples
//!
//! ```
//! use spacerocks::data::{constants, observatories};
//!
//! // Access physical constants
//! let c = constants::SPEED_OF_LIGHT;  // AU/day
//! let g = constants::GRAVITATIONAL_CONSTANT;
//!
//! // Get observatory coordinates
//! if let Some(&coords) = observatories::OBSERVATORIES.get("568") {
//!     let (lon, lat_cos, height_cos) = coords;  // Mauna Kea
//!     // lon ≈ 3.57 rad (-155.5°)
//!     // lat_cos and height_cos are direction cosines
//! }
//! ```


pub mod constants;
    pub use self::constants::*;

pub mod observatories;
    pub use self::observatories::*;