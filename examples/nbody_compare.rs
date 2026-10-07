//! Integrate a system read from stdin with one of the nbody integrators, for comparisons with
//! other codes (see `validation/rebound/compare.py`).
//!
//! Input: a header line `integrator dt n_steps n_out` (integrator: whfast | trace | ias15 |
//! ias15i | leapfrog; IAS15 integrates to `n_steps * dt` with its own adaptive steps, landing on
//! each output with `integrate`, or with `integrate_or_interpolate` for ias15i), then one particle
//! per line: `name mass x y z vx vy vz` (Msun, AU, AU/day).
//! Output: `energy t E` lines at `n_out` evenly spaced times, `state name x y z vx vy vz`
//! for every particle at the end, and `seconds s` for the time spent stepping.

use std::io::Read;
use std::time::Instant;

use spacerocks::nbody::{Integrator, Leapfrog, Simulation, Trace, WisdomHolman, IAS15};
use spacerocks::{SpaceRock, Time};

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let mut lines = input.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let dt: f64 = header[1].parse().unwrap();
    let n_steps: usize = header[2].parse().unwrap();
    let n_out: usize = header[3].parse().unwrap();
    let integrator: Box<dyn Integrator + Send + Sync> = match header[0] {
        "whfast" => Box::new(WisdomHolman::new(dt)),
        "trace" => Box::new(Trace::new(dt)),
        "ias15" | "ias15i" => Box::new(IAS15::new(dt)),
        "leapfrog" => Box::new(Leapfrog::new(dt)),
        other => panic!("unknown integrator {other}"),
    };

    let t0 = 2460000.5;
    let epoch = Time::new(t0, "tdb", "jd").unwrap();
    let mut sim = Simulation::new(&epoch, "J2000", "SSB").unwrap();
    sim.integrator = integrator;
    for line in lines {
        let f: Vec<&str> = line.split_whitespace().collect();
        let v: Vec<f64> = f[1..].iter().map(|x| x.parse().unwrap()).collect();
        let mut rock = SpaceRock::from_xyz(f[0], v[1], v[2], v[3], v[4], v[5], v[6], epoch.clone(), "J2000", "SSB").unwrap();
        rock.set_mass(v[0]);
        sim.add(rock).unwrap();
    }

    let per_out = (n_steps / n_out.max(1)).max(1);
    let mut seconds = 0.0;
    let mut taken = 0;
    println!("energy 0 {:.17e}", sim.energy());
    while taken < n_steps {
        let chunk = per_out.min(n_steps - taken);
        let start = Instant::now();
        if header[0] == "ias15" || header[0] == "ias15i" {
            // Adaptive: integrate to the output time instead of counting steps.
            let target = Time::new(t0 + (taken + chunk) as f64 * dt, "tdb", "jd").unwrap();
            if header[0] == "ias15i" {
                sim.integrate_or_interpolate(&target);
            } else {
                sim.integrate(&target);
            }
        } else {
            sim.steps(chunk);
        }
        seconds += start.elapsed().as_secs_f64();
        taken += chunk;
        if n_out > 0 {
            println!("energy {} {:.17e}", sim.epoch.tdb().jd() - t0, sim.energy());
        }
    }
    for (name, [x, y, z, vx, vy, vz]) in sim.names().iter().zip(sim.states()) {
        println!("state {} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}", name, x, y, z, vx, vy, vz);
    }
    println!("time {:.17e}", sim.epoch.tdb().jd() - t0);
    println!("seconds {seconds:.6}");
}
