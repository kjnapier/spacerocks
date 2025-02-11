//! N-body simulation management and control.
//! 
//! This module provides functionality for creating and managing n-body simulations,
//! including methods for:
//! - Setting up pre-configured solar system simulations
//! - Adding and removing particles
//! - Controlling integration
//! - Managing coordinate systems and reference frames
//! - Computing system properties like energy

pub mod simulation;

pub mod forces;

pub mod integrators;
    pub use self::integrators::Integrator;
    pub use self::integrators::Leapfrog;
    pub use self::integrators::IAS15;
    // pub use self::integrators::MVS;


pub use self::simulation::Simulation;
