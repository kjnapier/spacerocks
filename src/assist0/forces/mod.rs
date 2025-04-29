//! Force implementations for n-body simulations.
//!
//! This module provides various force models that can be used in n-body simulations:
//! - [`NewtonianGravity`]: Classical gravitational force between bodies
//! - [`SolarGR`]: Relativistic corrections for solar gravity
//! - [`SolarJ2`]: Perturbations from the Sun's oblateness
//!
//! Forces implement the [`Force`] trait which defines how they calculate accelerations
//! on a system of bodies.

pub mod force;
    pub use self::force::Force;

// pub mod drag;
//     pub use self::drag::Drag;

pub mod gravity;
    pub use self::gravity::NewtonianGravity;

// pub mod solar_gr;
//     pub use self::solar_gr::SolarGR;

// pub mod solar_j2;
//     pub use self::solar_j2::SolarJ2;