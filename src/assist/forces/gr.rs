//! Relativistic corrections, ported from ASSIST (`assist_additional_force_simple_GR`,
//! `assist_additional_force_potential_GR` and `assist_additional_force_eih_GR`), including the
//! terms for variational particles.

use nalgebra::{Matrix3, Vector3};

use crate::assist::constants::EphemerisConstants;
use crate::assist::forces::common::{add_jacobian, body, body_index, MAJOR_BODIES, SUN};
use crate::assist::forces::Force;
use crate::assist::SimulationState;

/// The Sun's post-Newtonian (PPN, beta = gamma = 1) one-body correction:
/// `a = GM/(c^2 r^3) [(4 GM/r - v^2) r + 4 (r.v) v]`, heliocentric position and velocity.
/// ASSIST's `GR_SIMPLE`.
#[derive(Debug, Clone, Copy)]
pub struct GrSimple {
    /// Speed of light squared, (AU/day)^2.
    pub c2: f64,
}

impl GrSimple {
    pub fn new(c: &EphemerisConstants) -> Self {
        GrSimple { c2: c.c_squared() }
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        let Some(is) = body_index(state, SUN) else { return };
        let (gm, xs, vs) = body(state, is);
        let c2 = self.c2;
        for (ip, p) in state.particles_1.iter_mut().enumerate() {
            let x = p.position - xs;
            let v = p.velocity - vs;
            let v2 = v.norm_squared();
            let r = x.norm();
            let a = 4.0 * gm / r - v2;
            let b = 4.0 * x.dot(&v);
            let prefac = gm / (r * r * r * c2);
            let f = a * x + b * v;
            p.acceleration += prefac * f;
            if !jac {
                continue;
            }
            let dpdr = -3.0 * prefac / r;
            let k = 4.0 * gm / (r * r);
            let mut dr = Matrix3::zeros();
            let mut dv = Matrix3::zeros();
            for i in 0..3 {
                for j in 0..3 {
                    let delta = if i == j { 1.0 } else { 0.0 };
                    dr[(i, j)] = dpdr * x[j] / r * f[i] + prefac * (delta * a - x[i] * (x[j] / r) * k + 4.0 * v[j] * v[i]);
                    dv[(i, j)] = prefac * (-2.0 * v[j] * x[i] + 4.0 * x[j] * v[i] + delta * b);
                }
            }
            add_jacobian(&mut state.partials[ip], &dr, Some(&dv));
        }
    }
}

impl Force for GrSimple {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}

/// A velocity-independent approximation to the Sun's relativistic correction,
/// `a = -6 (GM)^2 / (c^2 r^4) r_hat`. ASSIST's `GR_POTENTIAL`.
#[derive(Debug, Clone, Copy)]
pub struct GrPotential {
    pub c2: f64,
}

impl GrPotential {
    pub fn new(c: &EphemerisConstants) -> Self {
        GrPotential { c2: c.c_squared() }
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        let Some(is) = body_index(state, SUN) else { return };
        let (gm, xs, _) = body(state, is);
        for (ip, p) in state.particles_1.iter_mut().enumerate() {
            let x = p.position - xs;
            let r2 = x.norm_squared();
            let r = r2.sqrt();
            let prefac = -6.0 * gm * gm / (self.c2 * r2 * r2);
            p.acceleration += prefac * x;
            if !jac {
                continue;
            }
            let u = x / r;
            let mut dr = Matrix3::zeros();
            for i in 0..3 {
                for j in 0..3 {
                    dr[(i, j)] = if i == j { prefac } else { 0.0 } - 4.0 * prefac * u[i] * u[j];
                }
            }
            add_jacobian(&mut state.partials[ip], &dr, None);
        }
    }
}

impl Force for GrPotential {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}

/// Einstein–Infeld–Hoffmann (PPN, beta = gamma = 1) equations of motion for a test particle.
///
/// The relativistic terms are summed over the first `sources` of the Sun, Mercury, Venus,
/// Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune and Pluto (ASSIST's default is 1, the
/// Sun only); the potential sums inside each term run over all eleven. ASSIST's `GR_EIH`, its
/// default relativistic model.
#[derive(Debug, Clone, Copy)]
pub struct GrEih {
    pub c2: f64,
    pub sources: usize,
}

impl GrEih {
    pub fn new(c: &EphemerisConstants) -> Self {
        GrEih { c2: c.c_squared(), sources: 1 }
    }

    pub fn with_sources(mut self, sources: usize) -> Self {
        self.sources = sources.min(MAJOR_BODIES.len());
        self
    }

    fn apply(&self, state: &mut SimulationState, jac: bool) {
        // The eleven major bodies present in the simulation; the sources are those among the
        // first `sources` of them.
        let mut bodies: Vec<(f64, Vector3<f64>, Vector3<f64>)> = Vec::with_capacity(MAJOR_BODIES.len());
        let mut sources: Vec<usize> = Vec::new();
        for (n, &code) in MAJOR_BODIES.iter().enumerate() {
            if let Some(i) = body_index(state, code) {
                if n < self.sources {
                    sources.push(bodies.len());
                }
                bodies.push(body(state, i));
            }
        }
        let over_c2 = 1.0 / self.c2;
        let (beta, gamma) = (1.0f64, 1.0f64);

        // For each source j: sum_k GM_k / r_jk and its Newtonian acceleration from the other
        // major bodies. These do not depend on the particle.
        let per_source: Vec<(f64, Vector3<f64>)> = sources
            .iter()
            .map(|&j| {
                let xj = bodies[j].1;
                let mut term1 = 0.0;
                let mut aj = Vector3::zeros();
                for (k, &(gmk, xk, _)) in bodies.iter().enumerate() {
                    if k != j {
                        let djk = xj - xk;
                        let rjk2 = djk.norm_squared();
                        let rjk = rjk2.sqrt();
                        term1 += gmk / rjk;
                        aj -= gmk / (rjk2 * rjk) * djk;
                    }
                }
                (term1 * (-(2.0 * beta - 1.0) * over_c2), aj)
            })
            .collect();

        for (ip, p) in state.particles_1.iter_mut().enumerate() {
            let xi = p.position;
            let vi = p.velocity;
            let vi2 = vi.norm_squared();

            let mut acc = Vector3::zeros();
            let mut term7_sum = Vector3::zeros();
            let mut term8_sum = Vector3::zeros();
            let mut dr = Matrix3::zeros();
            let mut dv = Matrix3::zeros();
            let mut dterm7_dr = Matrix3::zeros();
            let mut dterm7_dv = Matrix3::zeros();
            let mut dterm8_dr = Matrix3::zeros();

            for (&j, &(term1, aj)) in sources.iter().zip(&per_source) {
                let (gmj, xj, vj) = bodies[j];
                let dij = xi - xj;
                let rij2 = dij.norm_squared();
                let rij = rij2.sqrt();
                let prefacij = gmj / (rij2 * rij);

                let term2 = gamma * over_c2 * vi2;
                let term3 = (1.0 + gamma) * over_c2 * vj.norm_squared();
                let vidotvj = vi.dot(&vj);
                let term4 = -2.0 * (1.0 + gamma) * over_c2 * vidotvj;
                let rijdotvj = dij.dot(&vj);
                let term5 = -1.5 * over_c2 * (rijdotvj * rijdotvj) / (rij * rij);

                let fvec = (2.0 + 2.0 * gamma) * vi - (1.0 + 2.0 * gamma) * vj;
                let f = dij.dot(&fvec);
                let vrel = vi - vj;
                term7_sum += prefacij * f * vrel;

                let mut term0 = 0.0;
                let mut dterm0 = Vector3::zeros();
                for &(gmk, xk, _) in bodies.iter() {
                    let dik = xi - xk;
                    let rik = dik.norm();
                    term0 += gmk / rik;
                    if jac {
                        dterm0 -= gmk / (rik * rik * rik) * dik;
                    }
                }
                term0 *= -2.0 * (beta + gamma) * over_c2;
                let term6 = -0.5 * over_c2 * dij.dot(&aj);
                let term8_fac = gmj / rij * (3.0 + 4.0 * gamma) / 2.0;
                term8_sum += term8_fac * aj;

                let factor = term0 + term1 + term2 + term3 + term4 + term5 + term6;
                acc += -prefacij * factor * dij;

                if !jac {
                    continue;
                }
                let dprefac = -3.0 * gmj / (rij2 * rij2 * rij) * dij; // d(prefacij)/dx_i
                let dterm0 = dterm0 * (-2.0 * (beta + gamma) * over_c2);
                let dterm2_dv = 2.0 * gamma * over_c2 * vi;
                let dterm4_dv = -2.0 * (1.0 + gamma) * over_c2 * vj;
                let term5_fac = 3.0 * over_c2 * rijdotvj / rij;
                let dterm5 = -term5_fac * (vj / rij - rijdotvj * dij / (rij * rij * rij));
                let dterm6 = -0.5 * over_c2 * aj;
                let dfactor_dr = dterm0 + dterm5 + dterm6;
                let dfactor_dv = dterm2_dv + dterm4_dv;

                // d/dx of (-prefacij * dij * factor)
                for a in 0..3 {
                    for b in 0..3 {
                        let delta = if a == b { 1.0 } else { 0.0 };
                        dr[(a, b)] += -dprefac[b] * dij[a] * factor - prefacij * factor * delta - prefacij * dij[a] * dfactor_dr[b];
                        dv[(a, b)] += -prefacij * dij[a] * dfactor_dv[b];
                    }
                }
                // term 7: prefacij * f * (vi - vj)
                let dfdv = (2.0 + 2.0 * gamma) * dij;
                for a in 0..3 {
                    for b in 0..3 {
                        let delta = if a == b { 1.0 } else { 0.0 };
                        dterm7_dr[(a, b)] += dprefac[b] * f * vrel[a] + prefacij * fvec[b] * vrel[a];
                        dterm7_dv[(a, b)] += prefacij * dfdv[b] * vrel[a] + prefacij * f * delta;
                    }
                }
                // term 8: gmj/rij (3+4 gamma)/2 * aj
                for a in 0..3 {
                    for b in 0..3 {
                        dterm8_dr[(a, b)] += -gmj * aj[a] / (rij * rij * rij) * dij[b] * (3.0 + 4.0 * gamma) / 2.0;
                    }
                }
            }
            p.acceleration += acc + term7_sum * over_c2 + term8_sum * over_c2;
            if jac {
                let dr = dr + (dterm7_dr + dterm8_dr) * over_c2;
                let dv = dv + dterm7_dv * over_c2;
                add_jacobian(&mut state.partials[ip], &dr, Some(&dv));
            }
        }
    }
}

impl Force for GrEih {
    fn apply_acceleration(&self, state: &mut SimulationState) {
        self.apply(state, false);
    }
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState) {
        self.apply(state, true);
    }
}
