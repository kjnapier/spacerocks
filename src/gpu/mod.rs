//! GPU integrators for large populations of test particles, with wgpu compute shaders (WGSL).
//!
//! Enabled with the `gpu` cargo feature. The integrators run the massive bodies on the CPU and
//! every test particle in its own GPU thread:
//!
//! - [`GpuWisdomHolman`]: the massive bodies take a block of Wisdom–Holman steps alone,
//!   recording what the test particles need from each step, and one kernel takes every test
//!   particle through the block (the same split as [`WisdomHolman`](crate::nbody::WisdomHolman)'s
//!   blocked stepping).
//! - [`GpuIAS15`]: the massive bodies are integrated with IAS15, keeping each step's
//!   Gauss–Radau polynomial, and each test particle runs its own IAS15, with its own adaptive
//!   timestep, through the bodies' trajectory.
//!
//! See [`Precision`] for the arithmetic: `DF64` works on every GPU, `F64` only where the
//! adapter has native f64, `F32` is not accurate enough to integrate with.

pub mod context;
pub use self::context::{GpuContext, GpuError, Precision};

pub mod wisdom_holman;
pub use self::wisdom_holman::GpuWisdomHolman;
