//! How far two-body motion from an N-body reference state drifts from N-body motion, as seen
//! from the Earth, for real MPCORB orbits (sets the checker's coarse-gate margin).
//! cargo run --release --example checker_twobody_error -- mpcorb_extended.json.gz
use spacerocks::batch::{BatchOptions};
use spacerocks::checker::{mpc, Catalog};
use spacerocks::{Origin, SpiceKernel};
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s};
use nalgebra::Vector3;

fn kepler(p: Vector3<f64>, v: Vector3<f64>, mu: f64, dt: f64) -> Vector3<f64> {
    let r = p.norm();
    let vr = v.dot(&p) / r;
    let alpha = 2.0 / r - v.norm_squared() / mu;
    let chi = solve_for_universal_anomaly(r, vr, alpha, mu, dt, 1e-12, 200).unwrap();
    let z = alpha * chi * chi;
    let f = 1.0 - chi * chi / r * stumpff_c(z);
    let g = dt - chi.powi(3) / mu.sqrt() * stumpff_s(z);
    p * f + v * g
}

fn main() {
    let path = std::env::args().nth(1).expect("path");
    let dir = std::path::PathBuf::from(std::env::var("SPACEROCKS_KERNELS").unwrap());
    let mut k = SpiceKernel::new();
    for f in ["de440s.bsp", "sb441-n16.bsp", "latest_leapseconds.tls"] { k.load(dir.join(f)).unwrap(); }
    let t = std::time::Instant::now();
    let els = mpc::read_mpcorb(std::path::Path::new(&path)).unwrap();
    println!("read {} orbits in {:.1} s", els.len(), t.elapsed().as_secs_f64());
    let t = std::time::Instant::now();
    let cat = Catalog::from_mpc_elements(&els, &k).unwrap();
    println!("catalog {} in {:.1} s", cat.len(), t.elapsed().as_secs_f64());
    // Sample: every 500th object plus all with q < 1.3.
    let mut idx: Vec<usize> = (0..cat.len()).step_by(500).collect();
    let neo: Vec<usize> = (0..cat.len()).filter(|&i| { let e = &els[i]; e.a * (1.0 - e.e) < 1.3 }).step_by(20).collect();
    idx.extend(&neo);
    let t_ref = 2460800.5 + 150.0; // snapshot epoch
    let dts = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0];
    let mut targets = vec![t_ref];
    targets.extend(dts.iter().map(|d| t_ref + d));
    let t = std::time::Instant::now();
    let s = cat.states_at(&idx, &targets, &k, &BatchOptions::default());
    println!("{} objects integrated in {:.1} s", idx.len(), t.elapsed().as_secs_f64());
    let mu = Origin::SUN.mu();
    let sun0 = k.state_au(10, 0, t_ref).unwrap();
    let m = targets.len();
    let classes = [("NEO q<1.3", 0.0, 1.3), ("q 1.3-5", 1.3, 5.0), ("q>5", 5.0, 1e9)];
    for (name, lo, hi) in classes {
        println!("{}", name);
        for (j, dt) in dts.iter().enumerate() {
            let sun = k.state_au(10, 0, targets[j + 1]).unwrap();
            let earth = k.state_au(399, 0, targets[j + 1]).unwrap();
            let mut errs = Vec::new();
            for (kk, &i) in idx.iter().enumerate() {
                let q = els[i].a * (1.0 - els[i].e);
                if q < lo || q >= hi { continue; }
                let r0 = s[kk * m];
                if !r0[0].is_finite() { continue; }
                let p = Vector3::new(r0[0] - sun0[0], r0[1] - sun0[1], r0[2] - sun0[2]);
                let v = Vector3::new(r0[3] - sun0[3], r0[4] - sun0[4], r0[5] - sun0[5]);
                let p2 = kepler(p, v, mu, *dt) + Vector3::new(sun[0], sun[1], sun[2]);
                let rn = s[kk * m + j + 1];
                let pn = Vector3::new(rn[0], rn[1], rn[2]);
                let e = Vector3::new(earth[0], earth[1], earth[2]);
                let a = (p2 - e).normalize();
                let b = (pn - e).normalize();
                errs.push((a - b).norm() * 206264.8);
            }
            errs.sort_by(|a, b| a.total_cmp(b));
            let q = |f: f64| errs[((errs.len() - 1) as f64 * f) as usize];
            println!("  dt {:5.2} d  n {:5}  median {:9.3}\"  99% {:9.3}\"  99.9% {:9.3}\"  max {:9.3}\"", dt, errs.len(), q(0.5), q(0.99), q(0.999), errs[errs.len() - 1]);
        }
    }
}
