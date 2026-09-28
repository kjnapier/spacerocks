//! Zonal harmonics of the Earth (J2, J3, J4) and the Sun (J2), ported from ASSIST
//! (`assist_additional_force_earth_J2J4` and `assist_additional_force_solar_J2`), including the
//! terms for variational particles.

use nalgebra::{Matrix3, Vector3};

use crate::assist::constants::EphemerisConstants;
use crate::assist::forces::common::{add_jacobian, body, body_index, pole_rotation, EARTH, SUN};
use crate::assist::forces::Force;
use crate::assist::SimulationState;

/// Earth's J2, J3 and J4 zonal harmonics.
///
/// As in ASSIST (and JPL Horizons), the Earth's pole is taken to be the J2000 pole, fixed in
/// time. Needs the Earth (NAIF 399) among the simulation's perturbers.
#[derive(Debug, Clone, Copy)]
pub struct EarthHarmonics {
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    /// Equatorial radius (AU).
    pub radius: f64,
    /// Pole right ascension and declination (radians).
    pub pole_ra: f64,
    pub pole_dec: f64,
}

impl EarthHarmonics {
    pub fn new(c: &EphemerisConstants) -> Self {
        EarthHarmonics {
            j2: c.j2e,
            j3: c.j3e,
            j4: c.j4e,
            radius: c.re_au(),
            pole_ra: 0.0,
            pole_dec: 90f64.to_radians(),
        }
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        let Some(ie) = body_index(state, EARTH) else { return };
        let (gm, xe, _) = body(state, ie);
        let rot = pole_rotation(self.pole_ra, self.pole_dec);
        let rot_t = rot.transpose();
        let re = self.radius;
        let (j2, j3, j4) = (self.j2, self.j3, self.j4);

        for p in state.particles_1.iter_mut() {
            let d0 = p.position - xe;
            let r2 = d0.norm_squared();
            let r = r2.sqrt();
            let d = rot * d0;
            let (dx, dy, dz) = (d.x, d.y, d.z);

            // J2
            let c2 = dz * dz / r2;
            let j2_prefac = 3.0 * j2 * re * re / r2 / r2 / r / 2.0;
            let j2_fac = 5.0 * c2 - 1.0;
            let mut res = Vector3::new(
                gm * j2_prefac * j2_fac * dx,
                gm * j2_prefac * j2_fac * dy,
                gm * j2_prefac * (j2_fac - 2.0) * dz,
            );
            // J3
            let j3_prefac = 5.0 * j3 * re * re * re / r2 / r2 / r / 2.0;
            let j3_fac = 3.0 - 7.0 * c2;
            res.x += -gm * j3_prefac * (1.0 / r2) * j3_fac * dx * dz;
            res.y += -gm * j3_prefac * (1.0 / r2) * j3_fac * dy * dz;
            res.z += -gm * j3_prefac * (6.0 * c2 - 7.0 * c2 * c2 - 0.6);
            // J4
            let j4_prefac = 5.0 * j4 * re * re * re * re / r2 / r2 / r2 / r / 8.0;
            let j4_fac = 63.0 * c2 * c2 - 42.0 * c2 + 3.0;
            res.x += gm * j4_prefac * j4_fac * dx;
            res.y += gm * j4_prefac * j4_fac * dy;
            res.z += gm * j4_prefac * (j4_fac + 12.0 - 28.0 * c2) * dz;

            p.acceleration += rot_t * res;

            if !jac {
                continue;
            }
            // J2
            let j2_fac2 = 7.0 * c2 - 1.0;
            let j2_fac3 = 35.0 * c2 * c2 - 30.0 * c2 + 3.0;
            let k2 = gm * j2_prefac;
            let dxdx = k2 * (j2_fac - 5.0 * j2_fac2 * dx * dx / r2);
            let dydy = k2 * (j2_fac - 5.0 * j2_fac2 * dy * dy / r2);
            let dzdz = k2 * (-1.0) * j2_fac3;
            let dxdy = k2 * (-5.0) * j2_fac2 * dx * dy / r2;
            let dydz = k2 * (-5.0) * (j2_fac2 - 2.0) * dy * dz / r2;
            let dxdz = k2 * (-5.0) * (j2_fac2 - 2.0) * dx * dz / r2;
            // J3
            let ct = dz / r;
            let j3_fac2 = 21.0 * (-3.0 * c2 + 1.0) / r2;
            let j3_fac3 = 3.0 * (-21.0 * c2 * c2 + 14.0 * c2 - 1.0) / r2;
            let j3_fac4 = (-63.0 * c2 * c2 + 70.0 * c2 - 15.0) * ct / r;
            let k3 = gm * j3_prefac;
            let dxdx3 = k3 * ct * (j3_fac2 * dx * dx - j3_fac) / r;
            let dydy3 = k3 * ct * (j3_fac2 * dy * dy - j3_fac) / r;
            let dzdz3 = k3 * j3_fac4;
            let dxdy3 = k3 * j3_fac2 * ct * dx * dy / r;
            let dydz3 = k3 * j3_fac3 * dy;
            let dxdz3 = k3 * j3_fac3 * dx;
            // J4
            let j4_fac2 = 33.0 * c2 * c2 - 18.0 * c2 + 1.0;
            let j4_fac3 = 33.0 * c2 * c2 - 30.0 * c2 + 5.0;
            let j4_fac4 = 231.0 * c2 * c2 * c2 - 315.0 * c2 * c2 + 105.0 * c2 - 5.0;
            let k4 = gm * j4_prefac;
            let dxdx4 = k4 * (j4_fac - 21.0 * j4_fac2 * dx * dx / r2);
            let dydy4 = k4 * (j4_fac - 21.0 * j4_fac2 * dy * dy / r2);
            let dzdz4 = k4 * (-3.0) * j4_fac4;
            let dxdy4 = k4 * (-21.0) * j4_fac2 * dx * dy / r2;
            let dydz4 = k4 * (-21.0) * j4_fac3 * dy * dz / r2;
            let dxdz4 = k4 * (-21.0) * j4_fac3 * dx * dz / r2;

            let (xx, yy, zz) = (dxdx + dxdx3 + dxdx4, dydy + dydy3 + dydy4, dzdz + dzdz3 + dzdz4);
            let (xy, yz, xz) = (dxdy + dxdy3 + dxdy4, dydz + dydz3 + dydz4, dxdz + dxdz3 + dxdz4);
            let m = Matrix3::new(xx, xy, xz, xy, yy, yz, xz, yz, zz);
            add_jacobian(p, &(rot_t * m * rot), None);
        }
    }
}

impl Force for EarthHarmonics {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}

/// The Sun's J2, with ASSIST's solar pole (RA 286.13°, Dec 63.87°).
#[derive(Debug, Clone, Copy)]
pub struct SolarJ2 {
    pub j2: f64,
    /// Solar radius (AU).
    pub radius: f64,
    pub pole_ra: f64,
    pub pole_dec: f64,
}

impl SolarJ2 {
    pub fn new(c: &EphemerisConstants) -> Self {
        SolarJ2 { j2: c.j2sun, radius: c.rsun_au(), pole_ra: 286.13f64.to_radians(), pole_dec: 63.87f64.to_radians() }
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        let Some(is) = body_index(state, SUN) else { return };
        let (gm, xs, _) = body(state, is);
        let rot = pole_rotation(self.pole_ra, self.pole_dec);
        let rot_t = rot.transpose();
        let rs = self.radius;

        for p in state.particles_1.iter_mut() {
            let d0 = p.position - xs;
            let r2 = d0.norm_squared();
            let r = r2.sqrt();
            let d = rot * d0;
            let (dx, dy, dz) = (d.x, d.y, d.z);

            let c2 = dz * dz / r2;
            let prefac = 3.0 * self.j2 * rs * rs / r2 / r2 / r / 2.0;
            let fac = 5.0 * c2 - 1.0;
            let fac2 = 7.0 * c2 - 1.0;
            let fac3 = 35.0 * c2 * c2 - 30.0 * c2 + 3.0;
            let res = Vector3::new(gm * prefac * fac * dx, gm * prefac * fac * dy, gm * prefac * (fac - 2.0) * dz);
            p.acceleration += rot_t * res;

            if !jac {
                continue;
            }
            let k = gm * prefac;
            let dxdx = k * (fac - 5.0 * fac2 * dx * dx / r2);
            let dydy = k * (fac - 5.0 * fac2 * dy * dy / r2);
            let dzdz = k * (-1.0) * fac3;
            let dxdy = k * (-5.0) * fac2 * dx * dy / r2;
            let dydz = k * (-5.0) * (fac2 - 2.0) * dy * dz / r2;
            let dxdz = k * (-5.0) * (fac2 - 2.0) * dx * dz / r2;
            let m = Matrix3::new(dxdx, dxdy, dxdz, dxdy, dydy, dydz, dxdz, dydz, dzdz);
            add_jacobian(p, &(rot_t * m * rot), None);
        }
    }
}

impl Force for SolarJ2 {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}
