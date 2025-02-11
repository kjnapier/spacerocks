//! Fundamental orbital mechanics data structures.
//!
//! This module provides the core structures used to represent orbits and states
//! in different coordinate systems:
//!
//! * [`KeplerOrbit`]: Represents orbits using Keplerian orbital elements
//! * [`StateVector`]: Represents position and velocity in Cartesian coordinates

pub mod keplerorbit;
pub use self::keplerorbit::KeplerOrbit;

pub mod statevector;
pub use self::statevector::StateVector;