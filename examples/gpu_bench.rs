//! Benchmarks and precision checks for the GPU integrators against the CPU ones, on the Sun and
//! giant planets with a population of Kuiper belt test particles (no kernels needed).
//!
//! Usage: `cargo run --release --features gpu --example gpu_bench -- <mode> [key=value ...]`
//!
//! - `wh`: Wisdom–Holman throughput. Keys: `n` (comma-separated particle counts, 1000,10000),
//!   `steps` (2000), `variants` (cpu,f32,df64,df64fast,f64,f64fast; `fast` = kicks in f32),
//!   `reps` (best of, 2). Prints test-particle steps per second and the largest relative
//!   position difference from the CPU after the run.
//! - `ias15`: IAS15 throughput. Keys: `n` (100,1000), `years` (10), `block` (365.25 days),
//!   `variants` (cpu,df64,f64), `pop` (`kbo`, or `centaurs`: planet crossers, whose encounters
//!   shorten the CPU's one shared timestep for everyone).
//! - `precision`: how far each variant drifts from the CPU over a long Wisdom–Holman run.
//!   Keys: `n` (200), `steps` (1e6), `outputs` (12, log-spaced), `planets` (1; 0 = the Sun
//!   alone, so orbits are regular and only round-off separates the variants), `variants`
//!   (f32,df64,df64fast,f64). Also runs the CPU from initial conditions moved by 1e-15
//!   (`cpu+1e-15`), the separation chaos alone produces. Prints, per output, the median and
//!   largest relative position difference and the median relative difference in semimajor axis.
//!
//! `SPACEROCKS_GPU=<name>` picks an adapter (e.g. `llvmpipe` for Mesa's software Vulkan).

#[path = "../tests/common/kuiper.rs"]
mod kuiper;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use kuiper::{kuiper, DT_WH};
use spacerocks::constants::GRAVITATIONAL_CONSTANT;
use spacerocks::gpu::{GpuContext, GpuIAS15, GpuWisdomHolman, Precision};
use spacerocks::nbody::{Integrator, Simulation, WisdomHolman, IAS15};
use spacerocks::Time;

fn list(opt: &HashMap<String, String>, key: &str, default: &str) -> Vec<String> {
    opt.get(key).map(|s| s.as_str()).unwrap_or(default).split(',').map(|s| s.trim().to_string()).collect()
}

fn num(opt: &HashMap<String, String>, key: &str, default: f64) -> f64 {
    opt.get(key).map(|s| s.parse().unwrap()).unwrap_or(default)
}

/// A variant's precision and whether its kicks are in f32, from names like `df64fast`.
fn variant(name: &str) -> Option<(Precision, bool)> {
    let fast = name.ends_with("fast");
    Precision::parse(name.trim_end_matches("fast")).map(|p| (p, fast))
}

fn wh_integrator(name: &str, ctx: &Arc<GpuContext>) -> Option<Box<dyn Integrator + Send + Sync>> {
    if name == "cpu" {
        return Some(Box::new(WisdomHolman::new(DT_WH)));
    }
    let (p, fast) = variant(name)?;
    match GpuWisdomHolman::new(DT_WH, ctx.clone(), p) {
        Ok(mut g) => {
            g.fast_kick = fast;
            g.min_work = 0;
            g.block = 8192;
            Some(Box::new(g))
        }
        Err(e) => {
            eprintln!("{name}: {e}");
            None
        }
    }
}

/// Relative position differences between the test particles of two simulations (sorted), and
/// the relative differences of their barycentric semimajor axes (sorted).
fn differences(a: &Simulation, b: &Simulation) -> (Vec<f64>, Vec<f64>) {
    let mu = GRAVITATIONAL_CONSTANT * a.masses().iter().sum::<f64>();
    let semi = |s: &[f64; 6]| {
        let r = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        let v2 = s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
        1.0 / (2.0 / r - v2 / mu)
    };
    let mut dx = vec![];
    let mut da = vec![];
    for ((sa, sb), &m) in a.states().iter().zip(b.states()).zip(a.masses()) {
        if m != 0.0 {
            continue;
        }
        let d = ((sa[0] - sb[0]).powi(2) + (sa[1] - sb[1]).powi(2) + (sa[2] - sb[2]).powi(2)).sqrt();
        let r = (sb[0] * sb[0] + sb[1] * sb[1] + sb[2] * sb[2]).sqrt();
        dx.push(d / r);
        da.push(((semi(sa) - semi(sb)) / semi(sb)).abs());
    }
    dx.sort_by(f64::total_cmp);
    da.sort_by(f64::total_cmp);
    (dx, da)
}

fn median(v: &[f64]) -> f64 {
    if v.is_empty() { f64::NAN } else { v[v.len() / 2] }
}

fn wh(opt: &HashMap<String, String>, ctx: &Arc<GpuContext>) {
    let steps = num(opt, "steps", 2000.0) as usize;
    let reps = num(opt, "reps", 2.0) as usize;
    println!("# Wisdom–Holman, Sun + 4 giants + n test particles, dt = {DT_WH:.2} d, {steps} steps, best of {reps}");
    println!("{:>8} {:>10} {:>12} {:>14} {:>12}", "n", "variant", "seconds", "steps/s", "max dx/x");
    for n in list(opt, "n", "1000,10000") {
        let n: usize = n.parse().unwrap();
        let mut reference: Option<Simulation> = None;
        for name in list(opt, "variants", "cpu,f32,df64,df64fast,f64,f64fast") {
            let Some(integrator) = wh_integrator(&name, ctx) else { continue };
            let base = kuiper(integrator, n, 1);
            let mut best = f64::INFINITY;
            let mut result = None;
            for _ in 0..reps {
                let mut sim = base.clone();
                let tic = Instant::now();
                sim.steps(steps);
                best = best.min(tic.elapsed().as_secs_f64());
                result = Some(sim);
            }
            let sim = result.unwrap();
            let dx = match &reference {
                Some(r) => differences(&sim, r).0.last().copied().unwrap_or(0.0),
                None => 0.0,
            };
            println!("{n:>8} {name:>10} {best:>12.4} {:>14.4e} {dx:>12.3e}", (n * steps) as f64 / best);
            if name == "cpu" {
                reference = Some(sim);
            }
        }
    }
}

fn ias15(opt: &HashMap<String, String>, ctx: &Arc<GpuContext>) {
    let years = num(opt, "years", 10.0);
    let block = num(opt, "block", 365.25);
    let pop = opt.get("pop").map(|s| s.as_str()).unwrap_or("kbo");
    println!("# IAS15, Sun + 4 giants + n test particles ({pop}), {years} years, GPU blocks of {block} days");
    println!("{:>8} {:>10} {:>12} {:>16} {:>12} {:>14} {:>10} {:>10}", "n", "variant", "seconds", "particle-yr/s", "max dx/x", "steps/particle", "cpu s", "gpu s");
    for n in list(opt, "n", "100,1000") {
        let n: usize = n.parse().unwrap();
        let mut reference: Option<Simulation> = None;
        for name in list(opt, "variants", "cpu,df64,f64") {
            let mut probe = None;
            let integrator: Box<dyn Integrator + Send + Sync> = if name == "cpu" {
                Box::new(IAS15::new(10.0))
            } else {
                let Some((p, _)) = variant(&name) else { continue };
                match GpuIAS15::new(block, ctx.clone(), p) {
                    Ok(g) => {
                        probe = Some(g.clone());
                        Box::new(g)
                    }
                    Err(e) => {
                        eprintln!("{name}: {e}");
                        continue;
                    }
                }
            };
            let mut sim = if pop == "centaurs" { kuiper::centaurs(integrator, n, 1) } else { kuiper(integrator, n, 1) };
            let end = Time::new(kuiper::T0 + years * 365.25, "tdb", "jd").unwrap();
            let tic = Instant::now();
            sim.integrate(&end);
            let secs = tic.elapsed().as_secs_f64();
            let dx = match &reference {
                Some(r) => differences(&sim, r).0.last().copied().unwrap_or(0.0),
                None => 0.0,
            };
            let (steps, cpu, gpu) = match &probe {
                Some(g) => {
                    let s = g.stats();
                    (format!("{:.1}", s.particle_steps / n as f64), format!("{:.3}", s.cpu_seconds), format!("{:.3}", s.gpu_wait_seconds))
                }
                None => ("-".into(), "-".into(), "-".into()),
            };
            println!("{n:>8} {name:>10} {secs:>12.3} {:>16.4e} {dx:>12.3e} {steps:>14} {cpu:>10} {gpu:>10}", n as f64 * years / secs);
            if name == "cpu" {
                reference = Some(sim);
            }
        }
    }
}

fn precision(opt: &HashMap<String, String>, ctx: &Arc<GpuContext>) {
    let n = num(opt, "n", 200.0) as usize;
    let steps = num(opt, "steps", 1e6) as usize;
    let outputs = num(opt, "outputs", 12.0) as usize;
    let planets = num(opt, "planets", 1.0) != 0.0;
    let make = |integrator: Box<dyn Integrator + Send + Sync>, shift: f64| kuiper::system(integrator, n, 1, planets, shift);
    let mut reference = make(Box::new(WisdomHolman::new(DT_WH)), 0.0);
    let mut runs: Vec<(String, Simulation)> = vec![("cpu+1e-15".to_string(), make(Box::new(WisdomHolman::new(DT_WH)), 1e-15))];
    for name in list(opt, "variants", "f32,df64,df64fast,f64") {
        if let Some(integrator) = wh_integrator(&name, ctx) {
            runs.push((name, make(integrator, 0.0)));
        }
    }
    println!("# precision: {n} test particles, {} planets, dt = {DT_WH:.2} d, up to {steps} steps ({:.3e} yr)", if planets { "Sun + 4 giant" } else { "Sun, no" }, steps as f64 * DT_WH / 365.25);
    println!("# relative position difference from the CPU (median, max), relative semimajor-axis difference (median)");
    println!("{:>12} {:>12} {:>10} {:>12} {:>12} {:>12}", "steps", "years", "variant", "med dx/x", "max dx/x", "med da/a");
    let mut done = 0;
    for k in 1..=outputs {
        let target = ((steps as f64).powf(k as f64 / outputs as f64)).round() as usize;
        if target <= done {
            continue;
        }
        let tic = Instant::now();
        reference.steps(target - done);
        let t_ref = tic.elapsed().as_secs_f64();
        for (name, sim) in runs.iter_mut() {
            let tic = Instant::now();
            sim.steps(target - done);
            let secs = tic.elapsed().as_secs_f64();
            let (dx, da) = differences(sim, &reference);
            println!("{target:>12} {:>12.4e} {name:>10} {:>12.3e} {:>12.3e} {:>12.3e}   ({secs:.1} s; cpu {t_ref:.1} s)", target as f64 * DT_WH / 365.25, median(&dx), dx.last().unwrap(), median(&da));
        }
        done = target;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args.first().cloned().unwrap_or_else(|| "wh".to_string());
    let opt: HashMap<String, String> = args.iter().skip(1).map(|a| {
        let (k, v) = a.split_once('=').unwrap_or_else(|| panic!("expected key=value, got {a}"));
        (k.to_string(), v.to_string())
    }).collect();
    let ctx = GpuContext::new().unwrap();
    println!("# GPU: {} (native f64: {}), CPU threads: {}", ctx.describe(), ctx.has_f64(), rayon::current_num_threads());
    match mode.as_str() {
        "wh" => wh(&opt, &ctx),
        "ias15" => ias15(&opt, &ctx),
        "precision" => precision(&opt, &ctx),
        _ => panic!("unknown mode {mode}: wh, ias15 or precision"),
    }
}
