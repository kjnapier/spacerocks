//! N-body simulation management and control.
//! 
//! This module provides functionality for creating and managing n-body simulations,
//! including methods for:
//! - Setting up pre-configured solar system simulations
//! - Adding and removing particles
//! - Controlling integration
//! - Managing coordinate systems and reference frames
//! - Computing system properties like energy

pub mod spice_simulation;
    pub use self::spice_simulation::SpiceSimulation;
    pub use self::spice_simulation::SimulationParticle;
    pub use self::spice_simulation::SimulationState;

    

pub mod constants;
    pub use self::constants::EphemerisConstants;

pub mod perturber_cache;
    pub use self::perturber_cache::PerturberCache;

pub mod forces;
pub mod integrators;
    pub use self::integrators::ias15::{AdaptiveMode, Summation, IAS15};
    pub use self::integrators::ias15::CoefficientSeptet;