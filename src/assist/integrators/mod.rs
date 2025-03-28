//! Numerical integrators for n-body simulations.
//!
//! This module provides different integration methods for advancing particles through time:
//! - [`IAS15`]: High-precision 15th order integrator with adaptive step size
//! - [`Leapfrog`]: Simple integrator with fixed step size
//!
//! Each integrator implements the [`Integrator`] trait which defines the core
//! functionality required for advancing particle states through time.

pub mod integrator;
    pub use self::integrator::Integrator;

pub mod ias15;
    pub use self::ias15::IAS15;

// pub mod mvs;
//     pub use self::mvs::MVS;