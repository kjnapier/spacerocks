//! Orbit determination and fitting.
//! 
//! The orbit module provides tools for orbit determination and fitting.
//! 
//! It contains implementations of:
//! - Gauss' method for initial orbit determination
//! - Levenberg-Marquardt optimization for orbit fitting
//! - Both analytic and numerical orbit propagation models
//! 
//! The module is structured around a few key traits and structs:
//! - `Model`: A trait for defining orbital models
//! - `LevenbergMarquardt`: Implementation of the LM optimization algorithm
//! - `AnalyticOrbitFitter` and `NumericalOrbitFitter`: Concrete orbit fitting implementations

pub mod gauss;
    pub use gauss::gauss;

pub mod fitter;
    pub use fitter::fit_orbit_lm;