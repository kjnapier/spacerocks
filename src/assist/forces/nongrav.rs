//! Non-gravitational accelerations (Marsden, Sekanina & Yeomans 1973), ported from ASSIST
//! (`assist_additional_force_non_gravitational`), including the terms for variational
//! particles and the partial derivatives with respect to A1, A2 and A3.

use nalgebra::Matrix3;

use crate::assist::forces::common::{add_jacobian, body, body_index, SUN};
use crate::assist::forces::Force;
use crate::assist::SimulationState;

/// Radial, transverse and normal accelerations `A_k g(r)` in the heliocentric RTN frame, with
/// `g(r) = alpha (r/r0)^-m (1 + (r/r0)^n)^-k`.
///
/// Each particle's (A1, A2, A3) come from its SpaceRock's `nongrav` property (AU/day^2);
/// particles without them feel no force. The default parameters (alpha = 1, r0 = 1 AU,
/// m = 2, n = 5.093, k = 0) give the inverse-square law JPL uses for asteroids; set them to the
/// water-ice sublimation values (alpha = 0.1112620426, r0 = 2.808, m = 2.15, n = 5.093,
/// k = 4.6142) for comets.
#[derive(Debug, Clone, Copy)]
pub struct NonGravitational {
    pub alpha: f64,
    pub r0: f64,
    pub m: f64,
    pub n: f64,
    pub k: f64,
}

impl Default for NonGravitational {
    fn default() -> Self {
        NonGravitational { alpha: 1.0, r0: 1.0, m: 2.0, n: 5.093, k: 0.0 }
    }
}

impl NonGravitational {
    /// The standard comet (water-ice sublimation) g(r) of Marsden et al. (1973).
    pub fn comet() -> Self {
        NonGravitational { alpha: 0.1112620426, r0: 2.808, m: 2.15, n: 5.093, k: 4.6142 }
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        if state.particles_1.iter().all(|p| p.nongrav == [0.0; 3]) {
            return;
        }
        let Some(is) = body_index(state, SUN) else { return };
        let (_, xs, vs) = body(state, is);
        let (alpha, r0, nm, nn, nk) = (self.alpha, self.r0, self.m, self.n, self.k);

        for p in state.particles_1.iter_mut() {
            let [a1, a2, a3] = p.nongrav;
            if a1 == 0.0 && a2 == 0.0 && a3 == 0.0 {
                continue;
            }
            let d = p.position - xs;
            let (dx, dy, dz) = (d.x, d.y, d.z);
            let r2 = d.norm_squared();
            let r = r2.sqrt();
            let g = alpha * (r / r0).powf(-nm) * (1.0 + (r / r0).powf(nn)).powf(-nk);

            let dv = p.velocity - vs;
            let (dvx, dvy, dvz) = (dv.x, dv.y, dv.z);
            let hx = dy * dvz - dz * dvy;
            let hy = dz * dvx - dx * dvz;
            let hz = dx * dvy - dy * dvx;
            let h = (hx * hx + hy * hy + hz * hz).sqrt();
            let tx = hy * dz - hz * dy;
            let ty = hz * dx - hx * dz;
            let tz = hx * dy - hy * dx;
            let t = (tx * tx + ty * ty + tz * tz).sqrt();

            p.acceleration.x += a1 * g * dx / r + a2 * g * tx / t + a3 * g * hx / h;
            p.acceleration.y += a1 * g * dy / r + a2 * g * ty / t + a3 * g * hy / h;
            p.acceleration.z += a1 * g * dz / r + a2 * g * tz / t + a3 * g * hz / h;

            if !jac {
                continue;
            }
            let r3 = r * r * r;
            let v2 = dvx * dvx + dvy * dvy + dvz * dvz;
            let rdotv = dx * dvx + dy * dvy + dz * dvz;
            let vdott = dvx * tx + dvy * ty + dvz * tz;

            let dgdr = (alpha / r0)
                * (-nm * (r / r0).powf(-nm - 1.0) * (1.0 + (r / r0).powf(nn)).powf(-nk)
                    + (r / r0).powf(-nm) * (-nk * nn) * (r / r0).powf(nn - 1.0) * (1.0 + (r / r0).powf(nn)).powf(-nk - 1.0));
            let dgx = dgdr * dx / r;
            let dgy = dgdr * dy / r;
            let dgz = dgdr * dz / r;

            let h3 = h * h * h;
            let (hxh3, hyh3, hzh3) = (hx / h3, hy / h3, hz / h3);
            let t3 = t * t * t;
            let (txt3, tyt3, tzt3) = (tx / t3, ty / t3, tz / t3);

            let dxdx = a1 * (dgx * dx / r + g * (1.0 / r - dx * dx / r3))
                + a2 * (dgx * tx / t + g * ((dx * dvx - rdotv) / t - txt3 * (2.0 * dx * vdott - rdotv * tx)))
                + a3 * (dgx * hx / h + g * (-hxh3) * (v2 * dx - rdotv * dvx));
            let dydy = a1 * (dgy * dy / r + g * (1.0 / r - dy * dy / r3))
                + a2 * (dgy * ty / t + g * ((dy * dvy - rdotv) / t - tyt3 * (2.0 * dy * vdott - rdotv * ty)))
                + a3 * (dgy * hy / h + g * (-hyh3) * (v2 * dy - rdotv * dvy));
            let dzdz = a1 * (dgz * dz / r + g * (1.0 / r - dz * dz / r3))
                + a2 * (dgz * tz / t + g * ((dz * dvz - rdotv) / t - tzt3 * (2.0 * dz * vdott - rdotv * tz)))
                + a3 * (dgz * hz / h + g * (-hzh3) * (v2 * dz - rdotv * dvz));
            let dxdy = a1 * (dgy * dx / r + g * (-dx * dy / r3))
                + a2 * (dgy * tx / t + g * ((2.0 * dy * dvx - dx * dvy) / t - txt3 * (2.0 * dy * vdott - rdotv * ty)))
                + a3 * (dgy * hx / h + g * (dvz / h - hxh3 * (v2 * dy - rdotv * dvy)));
            let dydx = a1 * (dgx * dy / r + g * (-dy * dx / r3))
                + a2 * (dgx * ty / t + g * ((2.0 * dx * dvy - dy * dvx) / t - tyt3 * (2.0 * dx * vdott - rdotv * tx)))
                + a3 * (dgx * hy / h + g * (-dvz / h - hyh3 * (v2 * dx - rdotv * dvx)));
            let dxdz = a1 * (dgz * dx / r + g * (-dx * dz / r3))
                + a2 * (dgz * tx / t + g * ((2.0 * dz * dvx - dx * dvz) / t - txt3 * (2.0 * dz * vdott - rdotv * tz)))
                + a3 * (dgz * hx / h + g * (-dvy / h - hxh3 * (v2 * dz - rdotv * dvz)));
            let dzdx = a1 * (dgx * dz / r + g * (-dz * dx / r3))
                + a2 * (dgx * tz / t + g * ((2.0 * dx * dvz - dz * dvx) / t - tzt3 * (2.0 * dx * vdott - rdotv * tx)))
                + a3 * (dgx * hz / h + g * (dvy / h - hzh3 * (v2 * dx - rdotv * dvx)));
            let dydz = a1 * (dgz * dy / r + g * (-dy * dz / r3))
                + a2 * (dgz * ty / t + g * ((2.0 * dz * dvy - dy * dvz) / t - tyt3 * (2.0 * dz * vdott - rdotv * tz)))
                + a3 * (dgz * hy / h + g * (dvx / h - hyh3 * (v2 * dz - rdotv * dvz)));
            let dzdy = a1 * (dgy * dz / r + g * (-dz * dy / r3))
                + a2 * (dgy * tz / t + g * ((2.0 * dy * dvz - dz * dvy) / t - tzt3 * (2.0 * dy * vdott - rdotv * ty)))
                + a3 * (dgy * hz / h + g * (-dvx / h - hzh3 * (v2 * dy - rdotv * dvy)));

            let dxdvx = a2 * g * ((dy * dy + dz * dz) / t - txt3 * r2 * tx) + a3 * g * (-hxh3 * (r2 * dvx - dx * rdotv));
            let dydvy = a2 * g * ((dx * dx + dz * dz) / t - tyt3 * r2 * ty) + a3 * g * (-hyh3 * (r2 * dvy - dy * rdotv));
            let dzdvz = a2 * g * ((dx * dx + dy * dy) / t - tzt3 * r2 * tz) + a3 * g * (-hzh3 * (r2 * dvz - dz * rdotv));
            let dxdvy = a2 * g * (-dy * dx / t - tyt3 * r2 * tx) + a3 * g * (-dz / h - hxh3 * (r2 * dvy - dy * rdotv));
            let dydvx = a2 * g * (-dx * dy / t - txt3 * r2 * ty) + a3 * g * (dz / h - hyh3 * (r2 * dvx - dx * rdotv));
            let dxdvz = a2 * g * (-dz * dx / t - tzt3 * r2 * tx) + a3 * g * (dy / h - hxh3 * (r2 * dvz - dz * rdotv));
            let dzdvx = a2 * g * (-dx * dz / t - txt3 * r2 * tz) + a3 * g * (-dy / h - hzh3 * (r2 * dvx - dx * rdotv));
            let dydvz = a2 * g * (-dz * dy / t - tzt3 * r2 * ty) + a3 * g * (-dx / h - hyh3 * (r2 * dvz - dz * rdotv));
            let dzdvy = a2 * g * (-dy * dz / t - tyt3 * r2 * tz) + a3 * g * (dx / h - hzh3 * (r2 * dvy - dy * rdotv));

            let dadr = Matrix3::new(dxdx, dxdy, dxdz, dydx, dydy, dydz, dzdx, dzdy, dzdz);
            let dadv = Matrix3::new(dxdvx, dxdvy, dxdvz, dydvx, dydvy, dydvz, dzdvx, dzdvy, dzdvz);
            add_jacobian(p, &dadr, Some(&dadv));

            // d(acceleration)/d(A1, A2, A3)
            let cols = [[g * dx / r, g * dy / r, g * dz / r], [g * tx / t, g * ty / t, g * tz / t], [g * hx / h, g * hy / h, g * hz / h]];
            for i in 0..3 {
                for k in 0..3 {
                    p.nongrav_partials[i][k] += cols[k][i];
                }
            }
        }
    }
}

impl Force for NonGravitational {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}
