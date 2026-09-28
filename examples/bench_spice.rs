//! Micro-benchmarks of common SPICE queries. Usage: cargo run --release --example bench_spice <kernel dir>
use rayon::prelude::*;
use spacerocks::spice::{et_from_jd, SpiceKernel};
use std::time::Instant;
fn main() {
    let d = std::env::args().nth(1).expect("usage: bench_spice <dir with de440s.bsp and earth_1962_240827_2124_combined.bpc>");
    let mut k = SpiceKernel::new();
    k.load(format!("{}/de440s.bsp", d)).unwrap();
    k.load(format!("{}/earth_1962_240827_2124_combined.bpc", d)).unwrap();
    let n = 200_000;
    let jds: Vec<f64> = (0..n).map(|i| 2460000.5 + (i as f64) * 0.013).collect();
    let ids = [10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9];
    let mut out = vec![[0.0; 6]; ids.len()];
    let mut acc = 0.0;
    let t = Instant::now(); for &jd in &jds { acc += k.state_au(399, 0, jd).unwrap()[0]; } let a = t.elapsed().as_nanos() as f64 / n as f64;
    let t = Instant::now(); for &jd in &jds { acc += k.state_au(301, 399, jd).unwrap()[0]; } let b = t.elapsed().as_nanos() as f64 / n as f64;
    let t = Instant::now(); for &jd in &jds[..n / 10] { k.barycentric_states_au(&ids, jd, &mut out).unwrap(); acc += out[3][0]; } let c = t.elapsed().as_nanos() as f64 / (n / 10) as f64;
    let t = Instant::now(); for &jd in &jds { acc += k.itrf93_to_j2000(et_from_jd(jd)).unwrap().rotation[0][0]; } let e = t.elapsed().as_nanos() as f64 / n as f64;
    let t = Instant::now(); let s: f64 = jds.par_iter().map(|&jd| k.state_au(399, 0, jd).unwrap()[0]).sum(); let f = t.elapsed().as_nanos() as f64 / n as f64;
    println!("NEW earth_wrt_ssb {:.0} ns | moon_wrt_earth {:.0} ns | 11-body barycentric {:.0} ns | earth pxform {:.0} ns | parallel earth_wrt_ssb {:.0} ns/call  ({})", a, b, c, e, f, acc + s);
}
