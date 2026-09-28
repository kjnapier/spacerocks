//! SPICE kernel support, implemented natively in Rust (no CSPICE dependency).
//!
//! * **SPK** (ephemerides): segment types 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21,
//!   with SPICE's precedence rules (later files and later segments win) and coverage-aware
//!   chaining between arbitrary targets and observers.
//! * **Binary PCK** (orientation): types 2, 3 and 20 (e.g. the high-precision Earth ITRF93 files).
//! * **Text kernels**: leapseconds, text PCKs (GM values, IAU rotation models), frame kernels
//!   (TK frames, body name definitions) and **meta-kernels** (`KERNELS_TO_LOAD`, path symbols).
//! * **Frames**: the 21 built-in inertial frames, built-in and kernel-defined PCK frames, and
//!   TK (fixed-offset) frames, with full state transformations (`sxform`) and rotations (`pxform`).
//!
//! Each [`SpiceKernel`] is an independent context; see the [`kernel`] module docs.
//!
//! Units: the SPICE-style methods ([`SpiceKernel::spkgeo`], [`SpiceKernel::state`],
//! [`SpiceKernel::sxform`], ...) use km, km/s and ephemeris time (TDB seconds past J2000).
//! The `*_au` helpers use AU, AU/day and TDB Julian dates like the rest of spacerocks.

mod builtin_bodies;
mod builtin_frames;
mod generic;
mod math;

pub mod bodies;
pub mod config;
pub mod daf;
pub mod error;
pub mod frames;
pub mod kernel;
pub mod pck;
pub mod spk;
pub mod text;

pub use bodies::SpiceBody;
pub use config::{Config as KernelConfig, KernelSpec};
pub use error::SpiceError;
pub use frames::{FrameInfo, StateTransform};
pub use kernel::{
    et_from_jd, iau_earth_rotation_delay, jd_from_et, km_to_au_state, Interval, KernelKind, KernelSummary, LoadedKernel, SegmentGroup, SpiceKernel, AU_KM, ITRF93_START_ET,
};
pub use pck::PckSegment;
pub use spk::SpkSegment;
pub use text::{KernelPool, PoolValue};
