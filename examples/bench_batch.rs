//! Accuracy and speed of the batch API (`spacerocks::batch`) against per-rock propagation.
//!
//!     SPACEROCKS_KERNELS=/path/to/kernels cargo run --release --example bench_batch
//!
//! The directory must hold `de440s.bsp`, `sb441-n16.bsp` and an Earth orientation file
//! (`earth_*.bpc`). Epochs are chosen inside the asteroid kernel's coverage.

use std::path::PathBuf;
use std::time::Instant;

use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use spacerocks::assist::{PerturberCache, SpiceSimulation};
use spacerocks::batch::{self, BatchOptions, Method};
use spacerocks::{Observatory, SpaceRock, SpiceKernel, Time};

const T0: f64 = 2460300.5;

fn kernel() -> SpiceKernel {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").expect("set SPACEROCKS_KERNELS"));
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).unwrap();
    k.load(dir.join("sb441-n16.bsp")).unwrap();
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if p.extension().map(|x| x == "bpc").unwrap_or(false) {
            k.load(&p).unwrap();
        }
    }
    k
}

#[derive(Clone, Copy, Debug)]
enum Class {
    Neo,
    MainBelt,
    Tno,
}

fn population(n: usize, classes: &[Class], seed: u64) -> Vec<(Class, SpaceRock)> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let t0 = Time::new(T0, "tdb", "jd").unwrap();
    (0..n)
        .map(|i| {
            let class = classes[i % classes.len()];
            let (q, e) = match class {
                Class::Neo => {
                    let e = rng.gen_range(0.2..0.6);
                    (rng.gen_range(0.8..1.3), e)
                }
                Class::MainBelt => {
                    let a: f64 = rng.gen_range(2.2..3.3);
                    let e = rng.gen_range(0.0..0.25);
                    (a * (1.0 - e), e)
                }
                Class::Tno => {
                    let a: f64 = rng.gen_range(35.0..50.0);
                    let e = rng.gen_range(0.0..0.3);
                    (a * (1.0 - e), e)
                }
            };
            let inc = rng.gen_range(0.0..0.5);
            let arg = rng.gen_range(0.0..std::f64::consts::TAU);
            let node = rng.gen_range(0.0..std::f64::consts::TAU);
            let f = rng.gen_range(0.0..std::f64::consts::TAU);
            let rock = SpaceRock::from_kepler(&format!("r{}", i), q, e, inc, arg, node, f, t0.clone(), "ECLIPJ2000", "SUN").unwrap();
            (class, rock)
        })
        .collect()
}

fn time<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let t = Instant::now();
    let r = f();
    (r, t.elapsed().as_secs_f64())
}

fn main() {
    let k = kernel();
    let threads = rayon::current_num_threads();
    println!("rayon threads: {}\n", threads);

    // --------------------------------------------------------------------------------------
    println!("== Perturber cache ==");
    let ids = SpiceSimulation::horizons_body_ids().unwrap();
    let (cache, dt) = time(|| PerturberCache::build(&k, &ids, T0 - 200.0, T0 + 400.0).unwrap());
    println!("build over 600 d: {:.1} ms", dt * 1e3);
    let mut iv = cache.intervals();
    iv.sort_by(|a, b| a.1.total_cmp(&b.1));
    println!("intervals (days): {:?}", iv.iter().map(|(id, l)| format!("{}:{:.3}", id, l)).collect::<Vec<_>>());
    let mut rng = rand::rngs::StdRng::seed_from_u64(1);
    let mut a = vec![[0.0; 6]; ids.len()];
    let mut b = vec![[0.0; 6]; ids.len()];
    let (mut ep, mut ev) = (0.0f64, 0.0f64);
    for _ in 0..20000 {
        let t = rng.gen_range(T0 - 200.0..T0 + 400.0);
        assert!(cache.states(&ids, t, &mut a));
        k.barycentric_states_au(&ids, t, &mut b).unwrap();
        for (x, y) in a.iter().zip(&b) {
            ep = ep.max(((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt());
            ev = ev.max(((x[3] - y[3]).powi(2) + (x[4] - y[4]).powi(2) + (x[5] - y[5]).powi(2)).sqrt());
        }
    }
    println!("max error at 20k random epochs: position {:.2e} AU ({:.2e} m), velocity {:.2e} AU/day", ep, ep * 1.495978707e11, ev);
    let n = 200_000;
    let (_, d_kernel) = time(|| {
        for i in 0..n {
            k.barycentric_states_au(&ids, T0 + i as f64 * 1e-3, &mut b).unwrap();
        }
    });
    let (_, d_cache) = time(|| {
        for i in 0..n {
            cache.states(&ids, T0 + i as f64 * 1e-3, &mut a);
        }
    });
    println!("27-body state lookup: kernel {:.2} us, cache {:.2} us", d_kernel / n as f64 * 1e6, d_cache / n as f64 * 1e6);
    let cache = std::sync::Arc::new(cache);
    let t0 = Time::new(T0, "tdb", "jd").unwrap();
    for np in [1usize, 16, 128] {
        let mut times = [0.0; 2];
        for (u, use_cache) in [false, true].into_iter().enumerate() {
            let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
            if use_cache {
                sim.set_perturber_cache(cache.clone());
            }
            for i in 0..np {
                let r = SpaceRock::from_kepler("x", 0.9 + i as f64 * 0.001, 0.4, 0.1, 1.0, 2.0, i as f64 * 0.05, t0.clone(), "J2000", "SSB").unwrap();
                sim.add(r).unwrap();
            }
            let (_, d) = time(|| {
                for j in 1..=12 {
                    sim.integrate_jd(T0 + 30.0 * j as f64, &k).unwrap();
                }
            });
            times[u] = d;
        }
        println!("{:>4} particles, 360 d: kernel {:.2} ms, cache {:.2} ms ({:.2}x)", np, times[0] * 1e3, times[1] * 1e3, times[0] / times[1]);
    }
    println!();

    // --------------------------------------------------------------------------------------
    println!("== Accuracy: batch propagate vs SpaceRock::propagate (+200 d) ==");
    let pop = population(300, &[Class::Neo, Class::MainBelt, Class::Tno], 7);
    let t1 = Time::new(T0 + 200.0, "tdb", "jd").unwrap();
    let reference: Vec<SpaceRock> = pop
        .par_iter()
        .map(|(_, r)| {
            let mut r = r.clone();
            r.propagate(&t1, &k).unwrap();
            r
        })
        .collect();
    for (label, opts) in [
        ("default (chunk 64)", BatchOptions::default()),
        ("with perturber cache", BatchOptions { perturber_cache: true, ..Default::default() }),
        ("chunk 1", BatchOptions { chunk_size: 1, ..Default::default() }),
    ] {
        let mut rocks: Vec<SpaceRock> = pop.iter().map(|(_, r)| r.clone()).collect();
        batch::propagate_batch(&mut rocks, &t1, &k, &opts).unwrap();
        let mut worst = [0.0f64; 3];
        for ((c, _), (a, b)) in pop.iter().zip(rocks.iter().zip(&reference)) {
            assert_eq!(a.origin, b.origin);
            assert_eq!(a.reference_plane, b.reference_plane);
            let d = (a.position - b.position).norm();
            let i = *c as usize;
            worst[i] = worst[i].max(d);
        }
        println!(
            "{:<22} max |dr| NEO {:.1e} AU, MBA {:.1e} AU, TNO {:.1e} AU  ({:.2} / {:.2} / {:.2} m)",
            label, worst[0], worst[1], worst[2], worst[0] * 1.496e11, worst[1] * 1.496e11, worst[2] * 1.496e11
        );
    }
    println!();

    // --------------------------------------------------------------------------------------
    println!("== Speed: propagate N rocks by +100 d ==");
    let t1 = Time::new(T0 + 100.0, "tdb", "jd").unwrap();
    for (label, classes, n) in [
        ("TNOs", vec![Class::Tno], 2000usize),
        ("main belt", vec![Class::MainBelt], 2000),
        ("mixed", vec![Class::Neo, Class::MainBelt, Class::Tno], 2000),
    ] {
        let pop: Vec<SpaceRock> = population(n, &classes, 11).into_iter().map(|(_, r)| r).collect();
        let (_, d_old) = time(|| {
            pop.clone().par_iter_mut().for_each(|r| r.propagate(&t1, &k).unwrap());
        });
        println!("{:<10} N={}  per-rock propagate (parallel): {:>8.1} ms  ({:.1} us/rock)", label, n, d_old * 1e3, d_old / n as f64 * 1e6);
        for chunk in [16usize, 64, 256, 1024] {
            let opts = BatchOptions { chunk_size: chunk, ..Default::default() };
            let mut rocks = pop.clone();
            let (_, d) = time(|| batch::propagate_batch(&mut rocks, &t1, &k, &opts).unwrap());
            println!(
                "           batch, chunk {:>4}: {:>8.1} ms  ({:.1} us/rock, {:.1}x)",
                chunk, d * 1e3, d / n as f64 * 1e6, d_old / d
            );
        }
    }
    println!();

    // --------------------------------------------------------------------------------------
    println!("== Ephemeris: N rocks x M epochs from W84 ==");
    let obs = Observatory::from_obscode("W84").unwrap();
    let m = 20;
    let epochs: Vec<Time> = (0..m).map(|j| Time::new(T0 - 60.0 + 6.0 * j as f64, "tdb", "jd").unwrap()).collect();
    let (observers, d_obs) = time(|| {
        epochs.par_iter().map(|t| obs.at(t, "J2000", "SSB", &k).unwrap()).collect::<Vec<_>>()
    });
    println!("{} observer states: {:.2} ms", m, d_obs * 1e3);

    let n = 1000;
    let pop: Vec<SpaceRock> = population(n, &[Class::MainBelt, Class::Tno], 3).into_iter().map(|(_, r)| r).collect();
    // Old way: for every epoch, propagate every rock from its epoch, then observe.
    let (old, d_old) = time(|| {
        pop.par_iter()
            .map(|r| {
                observers
                    .iter()
                    .map(|o| {
                        let mut x = r.clone();
                        x.propagate(&o.epoch, &k).unwrap();
                        x.change_reference_plane("J2000").unwrap();
                        x.to_ssb(&k).unwrap();
                        let ob = x.observe(o).unwrap();
                        (ob.ra(), ob.dec())
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    });
    println!("N={} M={}: propagate + observe per epoch: {:>8.1} ms", n, m, d_old * 1e3);
    for (label, opts) in [
        ("nbody", BatchOptions::default()),
        ("twobody", BatchOptions { method: Method::TwoBody, ..Default::default() }),
    ] {
        let (eph, d) = time(|| batch::ephemeris_rocks(&pop, &observers, &k, &opts).unwrap());
        let mut worst = 0.0f64;
        for i in 0..n {
            for j in 0..m {
                let a = eph.get(i, j);
                let (ra, dec) = old[i][j];
                let d = ((a.ra - ra) * dec.cos()).hypot(a.dec - dec);
                worst = worst.max(d);
            }
        }
        println!(
            "                 batch::ephemeris ({:<15}): {:>8.1} ms ({:.0}x), max offset from old {:.2e} mas",
            label, d * 1e3, d_old / d, worst.to_degrees() * 3.6e6
        );
    }
}
