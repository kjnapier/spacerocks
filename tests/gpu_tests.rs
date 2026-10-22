//! GPU integrators against the CPU ones. Runs on whatever adapter wgpu finds (in CI, Mesa's
//! software Vulkan driver); without one the tests print a note and pass.
#![cfg(feature = "gpu")]

mod common;

use std::sync::Arc;

use common::kuiper::{kuiper, DT_WH};
use spacerocks::gpu::{GpuContext, GpuWisdomHolman, Precision};
use spacerocks::nbody::{Simulation, WisdomHolman};

fn context() -> Option<Arc<GpuContext>> {
    match GpuContext::new() {
        Ok(ctx) => {
            eprintln!("GPU: {}", ctx.describe());
            Some(ctx)
        }
        Err(e) => {
            eprintln!("skipping, no GPU: {e}");
            None
        }
    }
}

fn precisions(ctx: &GpuContext) -> Vec<Precision> {
    let mut p = vec![];
    if ctx.has_f64() {
        p.push(Precision::F64);
    }
    p.push(Precision::DF64);
    p
}

/// Largest position and velocity differences between the two simulations' particles, relative
/// to each particle's distance and speed.
fn max_relative(a: &Simulation, b: &Simulation) -> (f64, f64) {
    a.particles().iter().zip(&b.particles()).fold((0.0, 0.0), |(dx, dv), (p, q)| {
        (f64::max(dx, (p.position - q.position).norm() / q.position.norm()), f64::max(dv, (p.velocity - q.velocity).norm() / q.velocity.norm()))
    })
}

#[test]
fn arithmetic_self_test() {
    let Some(ctx) = context() else { return };
    for p in [Precision::F32, Precision::DF64, Precision::F64] {
        if p == Precision::F64 && !ctx.has_f64() {
            continue;
        }
        let err = ctx.check(p).unwrap();
        eprintln!("{p:?}: worst relative error {err:e} ({:.2} units of round-off)", err / p.epsilon());
    }
}

#[test]
fn wisdom_holman_matches_cpu() {
    let Some(ctx) = context() else { return };
    for p in precisions(&ctx) {
        let mut gpu = GpuWisdomHolman::new(DT_WH, ctx.clone(), p).unwrap();
        gpu.min_work = 0;
        gpu.block = 300;
        let mut a = kuiper(Box::new(gpu), 200, 1);
        let mut b = kuiper(Box::new(WisdomHolman::new(DT_WH)), 200, 1);
        // Several blocks, the last one short.
        a.steps(1000);
        b.steps(1000);
        let (dx, dv) = max_relative(&a, &b);
        eprintln!("{p:?}: after 1000 steps, relative differences {dx:e} in position, {dv:e} in velocity");
        let tol = if p == Precision::F64 { 1e-12 } else { 1e-9 };
        assert!(dx < tol && dv < tol, "{p:?}: {dx:e} {dv:e}");
        // The massive bodies take exactly the CPU's steps.
        for (pa, pb) in a.particles().iter().zip(&b.particles()).take(5) {
            assert_eq!((pa.position, pa.velocity), (pb.position, pb.velocity), "{}", pa.name);
        }
        assert_eq!(a.epoch().tdb().jd(), b.epoch().tdb().jd());
    }
}

#[test]
fn wisdom_holman_small_calls_run_on_the_cpu() {
    let Some(ctx) = context() else { return };
    let gpu = GpuWisdomHolman::new(DT_WH, ctx, Precision::DF64).unwrap();
    let mut a = kuiper(Box::new(gpu), 10, 2);
    let mut b = kuiper(Box::new(WisdomHolman::new(DT_WH)), 10, 2);
    a.steps(100);
    b.steps(100);
    for (pa, pb) in a.particles().iter().zip(&b.particles()) {
        assert_eq!((pa.position, pa.velocity), (pb.position, pb.velocity), "{}", pa.name);
    }
}

#[test]
fn kepler_step_matches_cpu() {
    use spacerocks::constants::GRAVITATIONAL_CONSTANT;
    use spacerocks::transforms::universal_kepler_step;
    let Some(ctx) = context() else { return };
    let mu = GRAVITATIONAL_CONSTANT;
    let mut rnd = common::kuiper::Lcg(3);
    // Bound orbits from 1 to 50 au with e up to 0.9, short and long steps, and hyperbolic ones.
    let mut cases = vec![];
    for k in 0..2000 {
        let a: f64 = if k % 10 == 9 { -(1.0 + 5.0 * rnd.next()) } else { 1.0 + 49.0 * rnd.next() };
        let e: f64 = if a < 0.0 { 1.1 + 2.0 * rnd.next() } else { 0.9 * rnd.next() };
        let p = a * (1.0 - e * e);
        let nu: f64 = if e < 1.0 { std::f64::consts::TAU * rnd.next() } else { 0.9 * (2.0 * rnd.next() - 1.0) * (-1.0 / e).acos() };
        let r = p / (1.0 + e * nu.cos());
        let s = (mu / p).sqrt();
        let dt = [10.0, 216.6, 5000.0][k % 3] * (0.5 + rnd.next());
        cases.push(([r * nu.cos(), r * nu.sin(), 0.01 * r], [-s * nu.sin(), s * (e + nu.cos()), 0.0], mu, dt));
    }
    for p in precisions(&ctx).into_iter().chain([Precision::F32]) {
        let got = spacerocks::gpu::wisdom_holman::kepler_steps(&ctx, p, &cases);
        let mut worst: f64 = 0.0;
        let mut failed = 0;
        for ((r0, v0, gm, dt), (g, ok)) in cases.iter().zip(&got) {
            if !ok {
                failed += 1;
                continue;
            }
            let (r, v) = universal_kepler_step(&(*r0).into(), &(*v0).into(), *gm, *dt).unwrap();
            let er = ((g[0] - r.x).powi(2) + (g[1] - r.y).powi(2) + (g[2] - r.z).powi(2)).sqrt() / r.norm();
            let ev = ((g[3] - v.x).powi(2) + (g[4] - v.y).powi(2) + (g[5] - v.z).powi(2)).sqrt() / v.norm();
            worst = worst.max(er).max(ev);
        }
        eprintln!("{p:?}: Kepler step worst relative error {worst:e} ({:.0} units of round-off), {failed} of {} didn't converge", worst / p.epsilon(), cases.len());
        if p != Precision::F32 {
            assert!(worst < 1e4 * p.epsilon(), "{p:?}: {worst:e}");
            assert_eq!(failed, 0, "{p:?}");
        }
    }
}
