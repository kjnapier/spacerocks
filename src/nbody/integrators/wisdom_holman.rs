use crate::SpaceRock;
use crate::time::Time;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::integrators::{Integrator, Leapfrog};
use crate::nbody::forces::Force;
use crate::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s};

use nalgebra::Vector3;

/// A second-order Wisdom–Holman symplectic integrator in democratic heliocentric coordinates
/// (Duncan, Levison & Lee 1998).
///
/// The Hamiltonian is split into a Keplerian part about the central body, an interaction
/// part, and a "jump" part (from the central body's momentum), and each step is the
/// symmetric composition
///
/// kick(h/2) jump(h/2) kepler(h) jump(h/2) kick(h/2)
///
/// Positions are heliocentric (relative to the central body) and velocities barycentric.
/// The central body is the most massive particle. The Kepler drift is solved exactly with
/// universal variables, so for a test particle orbiting the central body alone the step is
/// exact, and with perturbers the energy error stays bounded instead of growing.
///
/// The kick uses the simulation's forces: the interaction acceleration on each body is its
/// total acceleration minus the central body's Keplerian pull, so extra forces (GR, J2, ...)
/// enter as perturbations. Velocity-dependent forces are evaluated at the start of each half
/// kick, which keeps the step explicit but not strictly symplectic.
///
/// The timestep is fixed. Close encounters between bodies other than the central one are not
/// handled, so choose a timestep well below the shortest orbital period (about 1/20 of it).
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct WisdomHolman {
    /// Current timestep in simulation time units
    pub timestep: f64,
}

impl WisdomHolman {
    /// Creates a new Wisdom–Holman integrator with the specified timestep.
    ///
    /// # Arguments
    ///
    /// * `timestep` - Fixed timestep to use for integration
    pub fn new(timestep: f64) -> WisdomHolman {
        WisdomHolman { timestep }
    }
}

/// Democratic heliocentric state: heliocentric positions `q`, barycentric velocities `u`,
/// and the barycenter's position and velocity. The central body's own `q` and `u` are unused.
#[derive(Clone)]
pub(crate) struct Democratic {
    pub(crate) central: usize,
    pub(crate) m_central: f64,
    pub(crate) m_total: f64,
    /// Masses, read once per step.
    pub(crate) m: Vec<f64>,
    pub(crate) q: Vec<Vector3<f64>>,
    pub(crate) u: Vec<Vector3<f64>>,
    pub(crate) x_cm: Vector3<f64>,
    pub(crate) v_cm: Vector3<f64>,
}

impl Democratic {
    pub(crate) fn from_particles(particles: &[SpaceRock], central: usize) -> Democratic {
        let m: Vec<f64> = particles.iter().map(|p| p.mass()).collect();
        let m_total: f64 = m.iter().sum();
        let mut x_cm = Vector3::zeros();
        let mut v_cm = Vector3::zeros();
        for (p, &mi) in particles.iter().zip(&m) {
            if mi != 0.0 {
                x_cm += mi * p.position;
                v_cm += mi * p.velocity;
            }
        }
        x_cm /= m_total;
        v_cm /= m_total;
        let xc = particles[central].position;
        let q = particles.iter().map(|p| p.position - xc).collect();
        let u = particles.iter().map(|p| p.velocity - v_cm).collect();
        Democratic { central, m_central: m[central], m_total, m, q, u, x_cm, v_cm }
    }

    /// Write inertial positions and velocities back into `particles`.
    pub(crate) fn to_particles(&self, particles: &mut [SpaceRock]) {
        let mut mq = Vector3::zeros();
        let mut mu = Vector3::zeros();
        for (i, &mi) in self.m.iter().enumerate() {
            if i != self.central && mi != 0.0 {
                mq += mi * self.q[i];
                mu += mi * self.u[i];
            }
        }
        let xc = self.x_cm - mq / self.m_total;
        let vc = self.v_cm - mu / self.m_central;
        for (i, p) in particles.iter_mut().enumerate() {
            if i == self.central {
                p.position = xc;
                p.velocity = vc;
            } else {
                p.position = self.q[i] + xc;
                p.velocity = self.u[i] + self.v_cm;
            }
        }
    }

    /// Kick the barycentric velocities by the forces minus the central body's Keplerian pull.
    /// The Newtonian pull between each pair in `skip` is left out too (TRACE moves those pairs
    /// into the drift), which assumes the forces include Newtonian gravity.
    pub(crate) fn kick(&mut self, particles: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>, h: f64, skip: &[(usize, usize)]) {
        self.to_particles(particles);
        let mut acc = vec![Vector3::zeros(); particles.len()];
        for force in forces {
            for (a, da) in acc.iter_mut().zip(force.calculate_acceleration(particles)) {
                *a += da;
            }
        }
        let gm = GRAVITATIONAL_CONSTANT * self.m_central;
        for i in 0..particles.len() {
            if i == self.central {
                continue;
            }
            let r = self.q[i].norm();
            let a_kepler = -gm * self.q[i] / (r * r * r);
            self.u[i] += h * (acc[i] - a_kepler);
        }
        for &(i, j) in skip {
            let d = self.q[j] - self.q[i];
            let r = d.norm();
            let g = GRAVITATIONAL_CONSTANT * d / (r * r * r);
            self.u[i] -= h * self.m[j] * g;
            self.u[j] += h * self.m[i] * g;
        }
    }

    pub(crate) fn jump(&mut self, h: f64) {
        let mut p = Vector3::zeros();
        for (i, &mi) in self.m.iter().enumerate() {
            if i != self.central && mi != 0.0 {
                p += mi * self.u[i];
            }
        }
        let dq = h * p / self.m_central;
        for i in 0..self.q.len() {
            if i != self.central {
                self.q[i] += dq;
            }
        }
    }

    pub(crate) fn kepler(&mut self, h: f64) {
        let gm = GRAVITATIONAL_CONSTANT * self.m_central;
        for i in 0..self.q.len() {
            if i != self.central {
                let (q, u) = kepler_drift(&self.q[i], &self.u[i], gm, h);
                self.q[i] = q;
                self.u[i] = u;
            }
        }
        self.x_cm += h * self.v_cm;
    }
}

/// Stumpff functions c0..c3 at `z`, by series at a reduced argument and the
/// quadruple-argument identities (no trigonometric calls).
#[inline]
fn stumpff_c0123(z: f64) -> [f64; 4] {
    let mut z = z;
    let mut n = 0;
    while z.abs() > 0.1 {
        z *= 0.25;
        n += 1;
    }
    // c2 = sum (-z)^j / (2j+2)!, c3 = sum (-z)^j / (2j+3)!, to z^6 (error < 1e-21 at |z| = 0.1)
    let c2 = 1.0 / 2.0 - z * (1.0 / 24.0 - z * (1.0 / 720.0 - z * (1.0 / 40320.0 - z * (1.0 / 3628800.0 - z * (1.0 / 479001600.0 - z / 87178291200.0)))));
    let c3 = 1.0 / 6.0 - z * (1.0 / 120.0 - z * (1.0 / 5040.0 - z * (1.0 / 362880.0 - z * (1.0 / 39916800.0 - z * (1.0 / 6227020800.0 - z / 1307674368000.0)))));
    let (mut c0, mut c1, mut c2, mut c3) = (1.0 - z * c2, 1.0 - z * c3, c2, c3);
    for _ in 0..n {
        c3 = (c2 + c0 * c3) * 0.25;
        c2 = 0.5 * c1 * c1;
        c1 *= c0;
        c0 = 2.0 * c0 * c0 - 1.0;
    }
    [c0, c1, c2, c3]
}

/// Advance a two-body orbit by `dt` with Gauss's f and g functions in universal variables
/// (Stumpff-function form, as in Wisdom & Hernandez 2015), solving Kepler's equation with
/// Halley's method. Falls back to the bracketing solver if that does not converge.
pub(crate) fn kepler_drift(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
    let r0n = r0.norm();
    let eta0 = r0.dot(v0);
    let beta = 2.0 * mu / r0n - v0.norm_squared();
    let zeta0 = mu - beta * r0n;

    // Whole periods of a bound orbit change nothing.
    let mut dt_red = dt;
    if beta > 0.0 {
        let period = std::f64::consts::TAU * mu / (beta * beta.sqrt());
        dt_red -= period * (dt / period).trunc();
    }

    // Initial guess: the short-step series, or the orbit-averaged rate 1/a for long steps.
    let mut x = if beta > 0.0 && dt_red.abs() * beta.sqrt() * beta / mu > 0.4 {
        dt_red * beta / mu
    } else {
        dt_red / r0n - 0.5 * eta0 * dt_red * dt_red / (r0n * r0n * r0n)
    };

    let mut converged = false;
    let mut last_dx = f64::INFINITY;
    for _ in 0..30 {
        let [c0, c1, c2, c3] = stumpff_c0123(beta * x * x);
        let (g1, g2, g3) = (x * c1, x * x * c2, x * x * x * c3);
        let f = r0n * x + eta0 * g2 + zeta0 * g3 - dt_red;
        let fp = r0n + eta0 * g1 + zeta0 * g2;
        let fpp = eta0 * c0 + zeta0 * g1;
        let dx = f / (fp - 0.5 * f * fpp / fp);
        x -= dx;
        if dx.abs() <= 2e-16 * x.abs() || (dx.abs() >= last_dx && dx.abs() <= 1e-12 * x.abs()) {
            converged = true;
            break;
        }
        last_dx = dx.abs();
    }
    if !converged || !x.is_finite() {
        return kepler_drift_bracketing(r0, v0, mu, dt);
    }

    let [_, c1, c2, c3] = stumpff_c0123(beta * x * x);
    let (g1, g2, g3) = (x * c1, x * x * c2, x * x * x * c3);
    let rn = r0n + eta0 * g1 + zeta0 * g2;
    let f = 1.0 - mu / r0n * g2;
    let g = dt_red - mu * g3;
    let fdot = -mu / (r0n * rn) * g1;
    let gdot = 1.0 - mu / rn * g2;
    (f * r0 + g * v0, fdot * r0 + gdot * v0)
}

/// Kepler drift with the robust bracketing universal-anomaly solver.
fn kepler_drift_bracketing(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
    let r0n = r0.norm();
    let vr0 = r0.dot(v0) / r0n;
    let alpha = 2.0 / r0n - v0.norm_squared() / mu;
    let sqrt_mu = mu.sqrt();
    let tol = 1e-15 * sqrt_mu * dt.abs();
    let chi = solve_for_universal_anomaly(r0n, vr0, alpha, mu, dt, tol, 100)
        .or_else(|_| solve_for_universal_anomaly(r0n, vr0, alpha, mu, dt, 1e3 * tol, 200))
        .expect("WisdomHolman: universal Kepler solve failed");

    let z = alpha * chi * chi;
    let (c, s) = (stumpff_c(z), stumpff_s(z));
    let f = 1.0 - chi * chi / r0n * c;
    let g = dt - chi * chi * chi / sqrt_mu * s;
    let r = f * r0 + g * v0;
    let rn = r.norm();
    let fdot = sqrt_mu / (rn * r0n) * chi * (z * s - 1.0);
    let gdot = 1.0 - chi * chi / rn * c;
    (r, fdot * r0 + gdot * v0)
}

/// The most massive particle (the first, on ties), if any has mass.
pub(crate) fn central_body(particles: &[SpaceRock]) -> Option<usize> {
    let best = particles.iter().enumerate().fold(None, |best: Option<(usize, f64)>, (i, p)| match best {
        Some((_, m)) if m >= p.mass() => best,
        _ => Some((i, p.mass())),
    });
    best.filter(|&(_, m)| m > 0.0).map(|(i, _)| i)
}

impl Integrator for WisdomHolman {
    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // The central body is the most massive particle (the first, on ties).
        let central = match central_body(particles) {
            Some(i) => i,
            // Nothing to orbit: fall back to a plain drift-kick-drift step.
            _ => return Leapfrog::new(self.timestep).step(particles, epoch, forces),
        };

        let h = self.timestep;
        let mut s = Democratic::from_particles(particles, central);
        s.kick(particles, forces, 0.5 * h, &[]);
        s.jump(0.5 * h);
        s.kepler(h);
        s.jump(0.5 * h);
        s.kick(particles, forces, 0.5 * h, &[]);
        s.to_particles(particles);

        *epoch += h;
        for particle in particles.iter_mut() {
            particle.epoch = epoch.clone();
        }
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kepler_drift_matches_bracketing_solver() {
        let mu = GRAVITATIONAL_CONSTANT;
        let cases = [
            (Vector3::new(1.0, 0.1, 0.0), Vector3::new(0.001, 0.017, 0.002), 30.0),   // near circular
            (Vector3::new(0.05, 0.0, 0.0), Vector3::new(0.0, 0.1, 0.0), 5.0),         // e ~ 0.98 at pericenter
            (Vector3::new(-4.9, 0.3, 0.0), Vector3::new(0.0, -0.0006, 0.0), 2000.0),  // near aphelion, long step
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.03, 0.0), 400.0),       // hyperbolic
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0172, 0.0), 3650.25),   // many periods
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0172, 0.0), -123.4),    // backwards
        ];
        for (r0, v0, dt) in cases {
            let (r, v) = kepler_drift(&r0, &v0, mu, dt);
            let (rb, vb) = kepler_drift_bracketing(&r0, &v0, mu, dt);
            assert!((r - rb).norm() < 1e-11 * rb.norm(), "r {r:?} vs {rb:?} ({r0:?}, {v0:?}, {dt})");
            assert!((v - vb).norm() < 1e-11 * vb.norm(), "v {v:?} vs {vb:?} ({r0:?}, {v0:?}, {dt})");
        }
    }
}
