//! Orbital mechanics transformation functions and solvers.
//!
//! This module provides functions for converting between different types of orbital
//! anomalies and solving Kepler's equation in universal variables. All angles are
//! in radians, distances in astronomical units (AU), and times in days.
//!
//! # Key Features
//!
//! - Conversion between mean, true, and conic anomalies for all orbit types
//! - Universal Kepler solver using Stumpff functions
//! - Light-time correction calculations
//!
//! # Anomaly Types
//!
//! The module handles three types of anomalies:
//! * Mean anomaly (M): Uniform angle increasing linearly with time
//! * True anomaly (ν): Actual angle of object from periapsis
//! * Conic anomaly: Generalized form that depends on orbit type:
//!   - For elliptical orbits (e < 1): Eccentric anomaly (E)
//!   - For hyperbolic orbits (e > 1): Hyperbolic anomaly (H)
//!   - For parabolic orbits (e = 1): Parabolic anomaly (B)
//!   - For circular orbits (e = 0): Equal to true anomaly
//!
//! # Units
//!
//! * Angles: radians
//! * Distances: astronomical units (AU)
//! * Time: days
//! * Velocities: AU/day
//!
//! # Examples
//!
//! ```
//! use spacerocks::transforms::{
//!     calc_conic_anomaly_from_mean_anomaly,
//!     calc_true_anomaly_from_conic_anomaly
//! };
//!
//! // Convert from mean to true anomaly for an elliptical orbit
//! let e = 0.7;  // eccentricity
//! let m = 0.8;  // mean anomaly in radians
//!
//! // First get conic (eccentric) anomaly
//! let e = calc_conic_anomaly_from_mean_anomaly(e, m).unwrap();
//!
//! // Then convert to true anomaly
//! let nu = calc_true_anomaly_from_conic_anomaly(e, e).unwrap();
//! ```

pub mod calc_conic_anomaly_from_true_anomaly;
    pub use self::calc_conic_anomaly_from_true_anomaly::calc_conic_anomaly_from_true_anomaly;

pub mod calc_conic_anomaly_from_mean_anomaly;
    pub use self::calc_conic_anomaly_from_mean_anomaly::calc_conic_anomaly_from_mean_anomaly;

pub mod correct_for_ltt;
    pub use self::correct_for_ltt::correct_for_ltt;

pub mod calc_mean_anomaly_from_conic_anomaly;
    pub use self::calc_mean_anomaly_from_conic_anomaly::calc_mean_anomaly_from_conic_anomaly;

pub mod calc_true_anomaly_from_conic_anomaly;
    pub use self::calc_true_anomaly_from_conic_anomaly::calc_true_anomaly_from_conic_anomaly;

pub mod calc_true_anomaly_from_mean_anomaly;
    pub use self::calc_true_anomaly_from_mean_anomaly::calc_true_anomaly_from_mean_anomaly;

pub mod stumpff;
    pub use self::stumpff::{stumpff_c, stumpff_s};

pub mod universal_kepler_solver;
    pub use self::universal_kepler_solver::solve_for_universal_anomaly;