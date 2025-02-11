//! Coordinate system types and transformations for dynamics calculations.
//!
//! This module provides the fundamental types for specifying coordinate systems:
//! - Origin points (Sun, Solar System Barycenter, or custom bodies)
//! - Reference planes (J2000, Ecliptic J2000, etc.)
//!
//! These components are essential for properly defining the position and motion of
//! celestial bodies in three-dimensional space.
//!
//! # Examples
//!
//! ```
//! use spacerocks::coordinates::{Origin, ReferencePlane};
//!
//! // Set up a heliocentric coordinate system
//! let origin = Origin::sun();
//! let plane = ReferencePlane::ECLIPJ2000;
//!
//! // Or create a custom system
//! let earth_centered = Origin::new_custom(0.000_000_000_889_954, "EARTH");
//! ```

pub mod origin;
    pub use self::origin::Origin;

pub mod reference_plane;
    pub use self::reference_plane::ReferencePlane;
