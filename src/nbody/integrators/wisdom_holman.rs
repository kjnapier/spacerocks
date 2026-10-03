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
struct Democratic {
    central: usize,
    m_central: f64,
    m_total: f64,
    q: Vec<Vector3<f64>>,
    u: Vec<Vector3<f64>>,
    x_cm: Vector3<f64>,
    v_cm: Vector3<f64>,
}

impl Democratic {
    fn from_particles(particles: &[SpaceRock], central: usize) -> Democratic {
        let m_total: f64 = particles.iter().map(|p| p.mass()).sum();
        let mut x_cm = Vector3::zeros();
        let mut v_cm = Vector3::zeros();
        for p in particles {
            x_cm += p.mass() * p.position;
            v_cm += p.mass() * p.velocity;
        }
        x_cm /= m_total;
        v_cm /= m_total;
        let xc = particles[central].position;
        let q = particles.iter().map(|p| p.position - xc).collect();
        let u = particles.iter().map(|p| p.velocity - v_cm).collect();
        Democratic { central, m_central: particles[central].mass(), m_total, q, u, x_cm, v_cm }
    }

    /// Write inertial positions and velocities back into `particles`.
    fn to_particles(&self, particles: &mut [SpaceRock]) {
        let mut mq = Vector3::zeros();
        let mut mu = Vector3::zeros();
        for (i, p) in particles.iter().enumerate() {
            if i != self.central {
                mq += p.mass() * self.q[i];
                mu += p.mass() * self.u[i];
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

    fn kick(&mut self, particles: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>, h: f64) {
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
    }

    fn jump(&mut self, particles: &[SpaceRock], h: f64) {
        let mut p = Vector3::zeros();
        for (i, rock) in particles.iter().enumerate() {
            if i != self.central {
                p += rock.mass() * self.u[i];
            }
        }
        let dq = h * p / self.m_central;
        for i in 0..self.q.len() {
            if i != self.central {
                self.q[i] += dq;
            }
        }
    }

    fn kepler(&mut self, h: f64) {
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

/// Advance a two-body orbit by `dt` with the universal-variable f and g functions.
fn kepler_drift(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
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

impl Integrator for WisdomHolman {
    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // The central body is the most massive particle (the first, on ties).
        let central = particles.iter().enumerate()
            .fold(None, |best: Option<(usize, f64)>, (i, p)| match best {
                Some((_, m)) if m >= p.mass() => best,
                _ => Some((i, p.mass())),
            });
        let central = match central {
            Some((i, m)) if m > 0.0 => i,
            // Nothing to orbit: fall back to a plain drift-kick-drift step.
            _ => return Leapfrog::new(self.timestep).step(particles, epoch, forces),
        };

        let h = self.timestep;
        let mut s = Democratic::from_particles(particles, central);
        s.kick(particles, forces, 0.5 * h);
        s.jump(particles, 0.5 * h);
        s.kepler(h);
        s.jump(particles, 0.5 * h);
        s.kick(particles, forces, 0.5 * h);
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
