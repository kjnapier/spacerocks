//! The Sun and the giant planets on their J2000 orbits (mean elements, so only roughly right)
//! with a population of Kuiper belt test particles, for the GPU tests and benchmarks. Needs no
//! kernels.

#![allow(dead_code)]

use spacerocks::constants::GRAVITATIONAL_CONSTANT;
use spacerocks::nbody::{Integrator, Simulation};
use spacerocks::{state, SpaceRock, Time};

pub const T0: f64 = 2451545.0;

/// Deterministic uniform numbers in [0, 1) (64-bit LCG).
pub struct Lcg(pub u64);

impl Lcg {
    pub fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// True anomaly from mean anomaly `m` (radians) on an ellipse of eccentricity `e`.
fn true_anomaly(m: f64, e: f64) -> f64 {
    let mut ea = m;
    for _ in 0..50 {
        ea -= (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
    }
    2.0 * ((1.0 + e).sqrt() * (0.5 * ea).sin()).atan2((1.0 - e).sqrt() * (0.5 * ea).cos())
}

fn rock(name: &str, mass: f64, s: [f64; 6]) -> SpaceRock {
    let mut r = SpaceRock::from_xyz(name, s[0], s[1], s[2], s[3], s[4], s[5], Time::new(T0, "tdb", "jd").unwrap(), "ECLIPJ2000", "SSB").unwrap();
    r.set_mass(mass);
    r
}

/// Sun (with the inner planets' mass) and the four giants, heliocentric, in the ecliptic.
pub fn giants(integrator: Box<dyn Integrator + Send + Sync>) -> Simulation {
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "ECLIPJ2000", "SSB").unwrap();
    sim.integrator = integrator;
    sim.add(rock("sun", 1.0 + 6.0e-6, [0.0; 6])).unwrap();
    // name, mass (solar), a, e, i, node, longitude of perihelion, mean longitude (deg)
    let planets = [
        ("jupiter", 9.547919e-4, 5.2026, 0.0485, 1.303, 100.46, 14.33, 34.35),
        ("saturn", 2.858860e-4, 9.5549, 0.0555, 2.489, 113.67, 93.06, 50.08),
        ("uranus", 4.366244e-5, 19.218, 0.0463, 0.773, 74.01, 173.01, 314.06),
        ("neptune", 5.151389e-5, 30.110, 0.0090, 1.770, 131.78, 48.12, 304.35),
    ];
    for (name, m, a, e, i, node, peri, l) in planets {
        let mu = GRAVITATIONAL_CONSTANT * (1.0 + m);
        let f = true_anomaly((l - peri as f64).to_radians(), e);
        let s = state::from_kepler(a * (1.0 - e), e, (i as f64).to_radians(), (peri - node as f64).to_radians(), (node as f64).to_radians(), f, mu).unwrap();
        sim.add(rock(name, m, s)).unwrap();
    }
    sim
}

/// [`giants`] plus `n` Kuiper belt test particles (a in 35..50 au, e < 0.2, i < 10 deg, random
/// angles, heliocentric elements), moved to the barycenter.
pub fn kuiper(integrator: Box<dyn Integrator + Send + Sync>, n: usize, seed: u64) -> Simulation {
    let mut sim = giants(integrator);
    let mu = GRAVITATIONAL_CONSTANT * (1.0 + 6.0e-6);
    let mut rnd = Lcg(seed);
    for k in 0..n {
        let a = 35.0 + 15.0 * rnd.next();
        let e = 0.2 * rnd.next();
        let i = (10.0 * rnd.next()).to_radians();
        let tau = std::f64::consts::TAU;
        let s = state::from_kepler(a * (1.0 - e), e, i, tau * rnd.next(), tau * rnd.next(), tau * rnd.next(), mu).unwrap();
        sim.add(rock(&format!("tp{k}"), 0.0, s)).unwrap();
    }
    sim.move_to_center_of_mass().unwrap();
    sim
}

/// Jupiter's period over 20, the usual Wisdom–Holman timestep for the giants (days).
pub const DT_WH: f64 = 4332.59 / 20.0;
