//! Stability of low-e, low-i test particles in the Kuiper belt under the giant planets, with
//! WHFast (`WisdomHolman`) and a long timestep.
//!
//! Usage: `kuiper_stability <de4xx.bsp> [n_particles] [a_min] [a_max] [years] [dt_days] [out_years]`
//!
//! The Sun (with the inner planets' masses added) and the four giant-planet barycenters come
//! from the kernel; the simulation is in the invariable plane about its barycenter. Test
//! particles are spaced evenly in a with e = 0.01 and i = 0.5 deg (barycentric elements about
//! the total mass) and random angles. Every `out_years` the barycentric a, e, i of each live
//! particle are printed as `el t_yr name a e inc_deg`. A particle is called unstable, printed as
//! `lost t_yr name a e q` and removed, once q < 30 au (Neptune-crossing) or e >= 1.

use std::time::Instant;

use spacerocks::constants::{GRAVITATIONAL_CONSTANT, MASSES};
use spacerocks::nbody::{Simulation, WisdomHolman};
use spacerocks::{state, SpaceRock, SpiceKernel, Time};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize, default: f64| args.get(i).map(|s| s.parse().unwrap()).unwrap_or(default);
    let kernel_path = &args[1];
    let n = arg(2, 21.0) as usize;
    let (a_min, a_max) = (arg(3, 35.0), arg(4, 39.0));
    let years = arg(5, 1.0e7);
    let dt = arg(6, 4332.59 / 20.0);
    let out_years = arg(7, 1000.0);

    let mut kernel = SpiceKernel::new();
    kernel.load(kernel_path).unwrap();
    let epoch = Time::new(2451545.0, "tdb", "jd").unwrap();
    let plane = "INVARIABLE";

    let mut sim = Simulation::new(&epoch, plane, "SSB").unwrap();
    sim.integrator = Box::new(WisdomHolman::new(dt));
    let inner: f64 = ["mercury barycenter", "venus barycenter", "earth barycenter", "mars barycenter"].iter().map(|b| MASSES[*b]).sum();
    let mut sun = SpaceRock::from_spice("sun", &epoch, plane, "SSB", &kernel).unwrap();
    sun.set_mass(MASSES["sun"] + inner);
    sim.add(sun).unwrap();
    for name in ["jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"] {
        sim.add(SpaceRock::from_spice(name, &epoch, plane, "SSB", &kernel).unwrap()).unwrap();
    }
    let mu = GRAVITATIONAL_CONSTANT * sim.masses().iter().sum::<f64>();

    // Deterministic angles (64-bit LCG), so runs are repeatable.
    let mut seed: u64 = 20261022;
    let mut uniform = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let (e0, i0) = (0.01, 0.5f64.to_radians());
    for k in 0..n {
        let a = if n == 1 { a_min } else { a_min + (a_max - a_min) * k as f64 / (n - 1) as f64 };
        let tau = std::f64::consts::TAU;
        let [x, y, z, vx, vy, vz] = state::from_kepler(a * (1.0 - e0), e0, i0, tau * uniform(), tau * uniform(), tau * uniform(), mu).unwrap();
        let rock = SpaceRock::from_xyz(&format!("tp{k:03}"), x, y, z, vx, vy, vz, epoch.clone(), plane, "SSB").unwrap();
        sim.add(rock).unwrap();
    }
    sim.move_to_center_of_mass().unwrap();

    let steps_per_out = ((out_years * 365.25 / dt).round() as usize).max(1);
    let n_out = (years / out_years).round() as usize;
    let e_start = sim.energy();
    println!("# n={n} a=[{a_min},{a_max}] years={years} dt_days={dt} out_years={out_years} steps={}", n_out * steps_per_out);

    let report = |sim: &mut Simulation, t_yr: f64| {
        let mut lost = Vec::new();
        for (name, (s, &m)) in sim.names().iter().zip(sim.states().iter().zip(sim.masses())) {
            if m != 0.0 {
                continue;
            }
            let e = state::eccentricity(s, mu);
            let a = state::semi_major_axis(s, mu);
            if e >= 1.0 || a * (1.0 - e) < 30.0 {
                println!("lost {t_yr:.1} {name} {a:.6} {e:.6} {:.6}", a * (1.0 - e));
                lost.push(name.clone());
            } else {
                println!("el {t_yr:.1} {name} {a:.6} {e:.6} {:.6}", state::inclination(s).to_degrees());
            }
        }
        for name in lost {
            sim.remove(&name).unwrap();
        }
    };

    report(&mut sim, 0.0);
    let start = Instant::now();
    for k in 1..=n_out {
        sim.steps(steps_per_out);
        let t_yr = sim.t() / 365.25;
        report(&mut sim, t_yr);
        if k % (n_out / 20).max(1) == 0 {
            let rel = ((sim.energy() - e_start) / e_start).abs();
            eprintln!("t = {:.3} Myr, {:.0} s elapsed, |dE/E| (incl. test-particle removals) = {rel:.2e}", t_yr / 1e6, start.elapsed().as_secs_f64());
        }
    }
    println!("# seconds {:.1}", start.elapsed().as_secs_f64());
}
