//! Identify detections with known objects, in the manner of the MPC's MPChecker.
//!
//! A [`Catalog`] holds orbits: the MPC's catalog of every minor planet (MPCORB, downloaded to
//! `~/.spacerocks/mpc`), orbits fitted with [`crate::orbfit`], or MPC `mpc_orb` records, each
//! with its uncertainty (a covariance, or for MPCORB the MPC's uncertainty parameter U).
//! [`check`] predicts every orbit at every detection and reports the objects each detection is
//! consistent with: those whose predicted position, with its uncertainty, lies within
//! [`CheckOptions::nsigma`] (Mahalanobis distance) of the detection, with its uncertainty.
//!
//! ```no_run
//! # use spacerocks::checker::{check, Catalog, CheckOptions};
//! # use spacerocks::orbfit::Astrometry;
//! # use spacerocks::SpiceKernel;
//! # fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
//! let kernel = SpiceKernel::defaults()?;
//! let mut catalog = Catalog::mpcorb(None, true, &kernel)?;      // ~1.4 million orbits
//! catalog.make_snapshot(2461000.5, &kernel, &Default::default(), 100_000);
//! catalog.save("mpcorb_2461000.srcat".as_ref())?;             // reuse: Catalog::load
//! # let detections = Astrometry::default();
//! let found = check(&catalog, &detections, &[], &kernel, &CheckOptions::default())?;
//! for m in found.matches.iter().filter(|m| m.consistent) {
//!     println!("detection {} is {} ({:.2} sigma)", m.detection, catalog.orbits.names[m.object], m.distance);
//! }
//! # Ok(()) }
//! ```

pub mod mpc;
    pub use mpc::{MpcElements, MpcOrb};

pub mod catalog;
    pub use catalog::Catalog;

pub mod check;
    pub use check::{check, CheckOptions, Checked, Match};
