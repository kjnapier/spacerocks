//! Stability of low-e, low-i test particles in the Kuiper belt under the giant planets, with
//! WHFast (`WisdomHolman`) and a long timestep, for runs long enough to need restarts.
//!
//! Usage: `kuiper_stability kernel=<de4xx.bsp> [key=value ...]`, with keys (defaults):
//! `n` (21), `a_min` (35), `a_max` (39), `e0` (0.01), `i0` (0.5, deg), `seed` (20261022),
//! `years` (1e7), `dt` (P_Jup/20 = 216.63 days), `sample` (1e4, years), `out` (1e3, years),
//! `run` (directory for `elements.txt` and `checkpoint.txt`; without it, elements go to stdout
//! and nothing is checkpointed), `checkpoint` (1e7, years).
//!
//! The Sun (with the inner planets' masses added) and the four giant-planet barycenters come
//! from the kernel; the simulation is in the invariable plane about its barycenter. Test
//! particles are spaced evenly in a, with barycentric e = e0 and i = i0 about the total mass and
//! random angles. The orbits are sampled every `sample` years. Every `out` years each live
//! particle gets a line `el t_yr name a e inc_deg a_min a_max e_max inc_max`: the elements now,
//! then the extremes over the samples since the last line. A particle is called unstable,
//! printed as `lost t_yr name a e q` and removed, at the first sample with q < 30 au or e >= 1.
//!
//! With `run`, the states are written to `checkpoint.txt` every `checkpoint` years, and a run
//! that finds a checkpoint carries on from it (the stepping is the same as an unbroken run's).
//! Raise `years` and rerun to extend a finished run.

use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use spacerocks::constants::{GRAVITATIONAL_CONSTANT, MASSES};
use spacerocks::nbody::{Simulation, WisdomHolman};
use spacerocks::{state, SpaceRock, SpiceKernel, State, Time};

const PLANE: &str = "INVARIABLE";
const J2000: f64 = 2451545.0;

/// Running extremes of one particle's elements between output lines.
#[derive(Clone, Copy)]
struct Extremes {
    a_min: f64,
    a_max: f64,
    e_max: f64,
    i_max: f64,
}

impl Extremes {
    fn new() -> Self {
        Extremes { a_min: f64::INFINITY, a_max: 0.0, e_max: 0.0, i_max: 0.0 }
    }

    fn add(&mut self, a: f64, e: f64, inc: f64) {
        self.a_min = self.a_min.min(a);
        self.a_max = self.a_max.max(a);
        self.e_max = self.e_max.max(e);
        self.i_max = self.i_max.max(inc);
    }
}

fn new_simulation(dt: f64) -> Simulation {
    let epoch = Time::new(J2000, "tdb", "jd").unwrap();
    let mut sim = Simulation::new(&epoch, PLANE, "SSB").unwrap();
    sim.integrator = Box::new(WisdomHolman::new(dt));
    sim
}

/// Giants from the kernel plus `n` test particles, about the barycenter.
fn initial(opt: &HashMap<String, String>, get: &dyn Fn(&str, f64) -> f64) -> Simulation {
    let mut kernel = SpiceKernel::new();
    kernel.load(&opt["kernel"]).unwrap();
    let epoch = Time::new(J2000, "tdb", "jd").unwrap();
    let mut sim = new_simulation(get("dt", 4332.59 / 20.0));
    let inner: f64 = ["mercury barycenter", "venus barycenter", "earth barycenter", "mars barycenter"].iter().map(|b| MASSES[*b]).sum();
    let mut sun = SpaceRock::from_spice("sun", &epoch, PLANE, "SSB", &kernel).unwrap();
    sun.set_mass(MASSES["sun"] + inner);
    sim.add(sun).unwrap();
    for name in ["jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"] {
        sim.add(SpaceRock::from_spice(name, &epoch, PLANE, "SSB", &kernel).unwrap()).unwrap();
    }
    let mu = GRAVITATIONAL_CONSTANT * sim.masses().iter().sum::<f64>();

    // Deterministic angles (64-bit LCG), so runs are repeatable.
    let mut seed = get("seed", 20261022.0) as u64;
    let mut uniform = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let n = get("n", 21.0) as usize;
    let (a_min, a_max) = (get("a_min", 35.0), get("a_max", 39.0));
    let (e0, i0) = (get("e0", 0.01), get("i0", 0.5).to_radians());
    let tau = std::f64::consts::TAU;
    for k in 0..n {
        let a = if n == 1 { a_min } else { a_min + (a_max - a_min) * k as f64 / (n - 1) as f64 };
        let [x, y, z, vx, vy, vz] = state::from_kepler(a * (1.0 - e0), e0, i0, tau * uniform(), tau * uniform(), tau * uniform(), mu).unwrap();
        sim.add(SpaceRock::from_xyz(&format!("tp{k:04}"), x, y, z, vx, vy, vz, epoch.clone(), PLANE, "SSB").unwrap()).unwrap();
    }
    sim.move_to_center_of_mass().unwrap();
    sim
}

/// Checkpoint: `t_days` and the starting energy on the first line, then `name mass x y z vx vy vz` as f64 bit patterns,
/// so a resumed run starts from exactly the same numbers.
fn write_checkpoint(path: &PathBuf, sim: &Simulation, t_days: f64, e_start: f64) {
    let tmp = path.with_extension("tmp");
    let mut f = File::create(&tmp).unwrap();
    writeln!(f, "{:016x} {:016x}", t_days.to_bits(), e_start.to_bits()).unwrap();
    for (name, (s, m)) in sim.names().iter().zip(sim.states().iter().zip(sim.masses())) {
        write!(f, "{name} {:016x}", m.to_bits()).unwrap();
        for v in s {
            write!(f, " {:016x}", v.to_bits()).unwrap();
        }
        writeln!(f).unwrap();
    }
    f.sync_all().unwrap();
    fs::rename(&tmp, path).unwrap();
}

fn read_checkpoint(path: &PathBuf, dt: f64) -> (Simulation, f64, f64) {
    let text = fs::read_to_string(path).unwrap();
    let bits = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).unwrap());
    let mut lines = text.lines();
    let head: Vec<f64> = lines.next().unwrap().split_whitespace().map(bits).collect();
    let epoch = Time::new(J2000, "tdb", "jd").unwrap();
    let mut sim = new_simulation(dt);
    for line in lines {
        // Names can hold spaces ("jupiter barycenter"); the last seven fields are numbers.
        let f: Vec<&str> = line.split_whitespace().collect();
        let (name, v) = f.split_at(f.len() - 7);
        let v: Vec<f64> = v.iter().map(|s| bits(s)).collect();
        let s: State = [v[1], v[2], v[3], v[4], v[5], v[6]];
        let mut rock = SpaceRock::from_xyz(&name.join(" "), s[0], s[1], s[2], s[3], s[4], s[5], epoch.clone(), PLANE, "SSB").unwrap();
        rock.set_mass(v[0]);
        sim.add(rock).unwrap();
    }
    (sim, head[0], head[1])
}

fn main() {
    let opt: HashMap<String, String> = std::env::args().skip(1).map(|a| {
        let (k, v) = a.split_once('=').unwrap_or_else(|| panic!("expected key=value, got {a}"));
        (k.to_string(), v.to_string())
    }).collect();
    let get = |k: &str, default: f64| opt.get(k).map(|s| s.parse().unwrap()).unwrap_or(default);
    let dt = get("dt", 4332.59 / 20.0);
    let years = get("years", 1.0e7);
    let steps_per_sample = ((get("sample", 1.0e4) * 365.25 / dt).round() as usize).max(1);
    let samples_per_out = ((get("out", 1.0e3) / get("sample", 1.0e4)).round() as usize).max(1);
    let samples_per_checkpoint = ((get("checkpoint", 1.0e7) / get("sample", 1.0e4)).round() as usize).max(1);
    let run = opt.get("run").map(PathBuf::from);

    // `t_days` counts days from J2000 across restarts; the simulation's own clock restarts at 0.
    let checkpoint = run.as_ref().map(|r| r.join("checkpoint.txt"));
    let (mut sim, t_days, e_start) = match &checkpoint {
        Some(path) if path.exists() => read_checkpoint(path, dt),
        _ => {
            let sim = initial(&opt, &get);
            let e = sim.energy();
            (sim, 0.0, e)
        }
    };
    let mut out: Box<dyn Write> = match &run {
        Some(r) => {
            fs::create_dir_all(r).unwrap();
            Box::new(OpenOptions::new().create(true).append(true).open(r.join("elements.txt")).unwrap())
        }
        None => Box::new(std::io::stdout()),
    };
    let mu = GRAVITATIONAL_CONSTANT * sim.masses().iter().sum::<f64>();
    let resumed = t_days > 0.0;
    if !resumed {
        writeln!(out, "# {}", std::env::args().skip(1).collect::<Vec<_>>().join(" ")).unwrap();
    }
    eprintln!("{} at {:.3} Myr with {} test particles", if resumed { "resuming" } else { "starting" }, t_days / 365.25e6,
              sim.masses().iter().filter(|&&m| m == 0.0).count());

    let mut extremes: HashMap<String, Extremes> = HashMap::new();
    let mut buffer = String::new();
    // Sample the elements; drop lost particles; with `line`, write `el` lines and reset extremes.
    let mut sample = |sim: &mut Simulation, t_yr: f64, line: bool, buffer: &mut String| {
        let mut lost = Vec::new();
        for (name, (s, &m)) in sim.names().iter().zip(sim.states().iter().zip(sim.masses())) {
            if m != 0.0 {
                continue;
            }
            let (a, e, inc) = (state::semi_major_axis(s, mu), state::eccentricity(s, mu), state::inclination(s).to_degrees());
            if e >= 1.0 || a * (1.0 - e) < 30.0 {
                buffer.push_str(&format!("lost {t_yr:.1} {name} {a:.6} {e:.6} {:.6}\n", a * (1.0 - e)));
                lost.push(name.clone());
                extremes.remove(name);
                continue;
            }
            let x = extremes.entry(name.clone()).or_insert_with(Extremes::new);
            x.add(a, e, inc);
            if line {
                buffer.push_str(&format!("el {t_yr:.1} {name} {a:.6} {e:.6} {inc:.5} {:.6} {:.6} {:.6} {:.5}\n", x.a_min, x.a_max, x.e_max, x.i_max));
                *x = Extremes::new();
            }
        }
        for name in lost {
            sim.remove(&name).unwrap();
        }
    };

    if !resumed {
        sample(&mut sim, 0.0, true, &mut buffer);
    }
    let start = Instant::now();
    let total_samples = (years * 365.25 / (steps_per_sample as f64 * dt)).round() as usize;
    let mut t_days = t_days;
    let mut k = (t_days / (steps_per_sample as f64 * dt)).round() as usize;
    while k < total_samples {
        sim.steps(steps_per_sample);
        k += 1;
        t_days = (k * steps_per_sample) as f64 * dt;
        sample(&mut sim, t_days / 365.25, k % samples_per_out == 0, &mut buffer);
        if k % samples_per_checkpoint == 0 || k == total_samples {
            out.write_all(buffer.as_bytes()).unwrap();
            out.flush().unwrap();
            buffer.clear();
            if let Some(path) = &checkpoint {
                write_checkpoint(path, &sim, t_days, e_start);
            }
            let left = sim.masses().iter().filter(|&&m| m == 0.0).count();
            eprintln!("t = {:.1} Myr, {:.0} s elapsed, {left} test particles left, |dE/E| = {:.2e}",
                      t_days / 365.25e6, start.elapsed().as_secs_f64(), ((sim.energy() - e_start) / e_start).abs());
        }
    }
}
