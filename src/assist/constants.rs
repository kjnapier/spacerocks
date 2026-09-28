//! Physical constants used by the force models.

use crate::spice::SpiceKernel;

/// Constants of the planetary ephemeris used by the non-Newtonian force models.
///
/// JPL planetary ephemerides record the constants they were integrated with in the SPK
/// comment area ("Initial conditions and constants used for integration"); ASSIST reads them
/// from there, and so does [`EphemerisConstants::from_kernel`]. Anything not found falls back
/// to the DE440 value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EphemerisConstants {
    /// Astronomical unit (km).
    pub au_km: f64,
    /// Speed of light (km/s).
    pub clight_km_s: f64,
    /// Earth zonal harmonics (unnormalized) and equatorial radius (km).
    pub j2e: f64,
    pub j3e: f64,
    pub j4e: f64,
    pub re_km: f64,
    /// Solar J2 and radius (km).
    pub j2sun: f64,
    pub asun_km: f64,
}

impl EphemerisConstants {
    /// DE440/DE441 values.
    pub const DE440: EphemerisConstants = EphemerisConstants {
        au_km: 1.495_978_707e8,
        clight_km_s: 2.997_924_58e5,
        j2e: 1.082_625_39e-3,
        j3e: -2.532_41e-6,
        j4e: -1.619_898_000_000_000_1e-6,
        re_km: 6.378_136_6e3,
        j2sun: 2.196_139_151_652_982_5e-7,
        asun_km: 6.96e5,
    };

    /// Read the constants from the comment areas of the loaded SPK files (highest priority
    /// first), falling back to DE440 for anything missing.
    pub fn from_kernel(kernel: &SpiceKernel) -> EphemerisConstants {
        let d = Self::DE440;
        let get = |k: &str, default: f64| kernel.integration_constant(k).unwrap_or(default);
        EphemerisConstants {
            au_km: get("AU", d.au_km),
            clight_km_s: get("CLIGHT", d.clight_km_s),
            j2e: get("J2E", d.j2e),
            j3e: get("J3E", d.j3e),
            j4e: get("J4E", d.j4e),
            re_km: get("RE", d.re_km),
            j2sun: get("J2SUN", d.j2sun),
            asun_km: get("ASUN", d.asun_km),
        }
    }

    /// Speed of light squared in (AU/day)^2.
    pub fn c_squared(&self) -> f64 {
        let c = self.clight_km_s / self.au_km * 86400.0;
        c * c
    }

    /// Earth's equatorial radius in AU.
    pub fn re_au(&self) -> f64 {
        self.re_km / self.au_km
    }

    /// Solar radius in AU.
    pub fn rsun_au(&self) -> f64 {
        self.asun_km / self.au_km
    }
}
