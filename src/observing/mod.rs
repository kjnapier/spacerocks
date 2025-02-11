//! Handling and manipulating astronomical observations.
//!
//! The observing module provides functionality for handling astronomical observations
//! and observers. It contains three main components:
//! 
//! - Observatory: Represents physical observation locations
//! - Observer: Represents an observatory at a specific time
//! - Observation: Handles different types of astronomical measurements
//!
//! The module supports both ground-based and space-based observations, as well as
//! various observation types including astrometry, radar, and streak observations.

pub mod observatory;
    pub use observatory::Observatory;

pub mod observer;
    pub use observer::Observer;

pub mod observation;
    pub use observation::{Observation};