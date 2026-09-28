//! Force models for [`crate::assist::SpiceSimulation`], following ASSIST (Holman et al. 2023).
//!
//! - [`NewtonianGravity`]: point-mass gravity of the SPICE perturbers
//! - [`EarthHarmonics`], [`SolarJ2`]: zonal harmonics of the Earth (J2, J3, J4) and the Sun (J2)
//! - [`GrEih`], [`GrSimple`], [`GrPotential`]: relativistic corrections
//! - [`NonGravitational`]: Marsden-model non-gravitational accelerations (A1, A2, A3)
//!
//! Every force also supplies the partial derivatives of its acceleration (with respect to
//! position, velocity and, for [`NonGravitational`], the A parameters) used by variational
//! particles. [`assist_default_forces`] returns ASSIST's default set.

pub mod force;
    pub use self::force::Force;

pub mod common;

pub mod gravity;
    pub use self::gravity::NewtonianGravity;

pub mod harmonics;
    pub use self::harmonics::{EarthHarmonics, SolarJ2};

pub mod gr;
    pub use self::gr::{GrEih, GrPotential, GrSimple};

pub mod nongrav;
    pub use self::nongrav::NonGravitational;

use crate::assist::constants::EphemerisConstants;

/// ASSIST's default force model, in ASSIST's order of application (smallest first):
/// non-gravitational forces (active only for particles with A1/A2/A3 set), Earth J2–J4,
/// solar J2, Einstein–Infeld–Hoffmann GR from the Sun, and Newtonian gravity.
pub fn assist_default_forces(c: &EphemerisConstants) -> Vec<Box<dyn Force + Send + Sync>> {
    vec![
        Box::new(NonGravitational::default()),
        Box::new(EarthHarmonics::new(c)),
        Box::new(SolarJ2::new(c)),
        Box::new(GrEih::new(c)),
        Box::new(NewtonianGravity),
    ]
}
