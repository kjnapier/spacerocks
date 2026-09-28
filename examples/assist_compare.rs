//! Driver for `validation/assist/compare.py`: integrates test particles with a chosen set of
//! forces and prints final states (and variational particles).
//!
//! stdin: a header line `t0 t1 forces epsilon [var]` (TDB JDs; forces a comma-separated list of
//! newton, ng, earth, sunj2, eih, eih11, grsimple, grpot), then one line per particle
//! `x y z vx vy vz A1 A2 A3` (barycentric J2000, AU, AU/day, AU/day^2).
//! Kernels: `de440s.bsp` and `sb441-n16.bsp` in `SPACEROCKS_KERNELS`.
use spacerocks::assist::forces::*;
use spacerocks::assist::{EphemerisConstants, SpiceSimulation, IAS15};
use spacerocks::{SpaceRock, SpiceKernel, Time};
use std::io::Read;
fn main() {
    let mut k = SpiceKernel::new();
    let dir = std::path::PathBuf::from(std::env::var("SPACEROCKS_KERNELS").expect("set SPACEROCKS_KERNELS"));
    k.load(dir.join("de440s.bsp")).unwrap();
    k.load(dir.join("sb441-n16.bsp")).unwrap();
    let c = EphemerisConstants::from_kernel(&k);
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let mut lines = input.lines();
    let head: Vec<String> = lines.next().unwrap().split_whitespace().map(String::from).collect();
    let (t0, t1, cfg, eps): (f64, f64, &str, f64) = (head[0].parse().unwrap(), head[1].parse().unwrap(), &head[2], head[3].parse().unwrap());
    let var = head.get(4).map(|s| s == "var").unwrap_or(false);
    for line in lines {
        let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
        if v.len() < 9 { continue; }
        let t = Time::new(t0, "tdb", "jd").unwrap();
        let mut r = SpaceRock::from_xyz("p", v[0], v[1], v[2], v[3], v[4], v[5], t.clone(), "J2000", "SSB").unwrap();
        r.set_nongrav(v[6], v[7], v[8]);
        let mut sim = SpiceSimulation::horizons(&t, &k).unwrap();
        let mut ias = IAS15::new(0.001); ias.epsilon = eps; sim.integrator = Box::new(ias);
        let mut forces: Vec<Box<dyn Force + Send + Sync>> = vec![];
        for f in cfg.split(',') {
            match f {
                "newton" => forces.push(Box::new(NewtonianGravity)),
                "ng" => forces.push(Box::new(NonGravitational::default())),
                "earth" => forces.push(Box::new(EarthHarmonics::new(&c))),
                "sunj2" => forces.push(Box::new(SolarJ2::new(&c))),
                "eih" => forces.push(Box::new(GrEih::new(&c))),
                "eih11" => forces.push(Box::new(GrEih::new(&c).with_sources(11))),
                "grsimple" => forces.push(Box::new(GrSimple::new(&c))),
                "grpot" => forces.push(Box::new(GrPotential::new(&c))),
                _ => panic!("{}", f),
            }
        }
        sim.set_forces(forces);
        sim.add(r).unwrap();
        if var { sim.add_full_variation("p").unwrap(); for a in ["A1","A2","A3"] { sim.add_variation(a, "p").unwrap(); } }
        sim.integrate_jd(t1, &k).unwrap();
        let p = &sim.state.particles[0];
        print!("{:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}", p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z);
        for vp in &sim.state.variational_particles { print!(" {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}", vp.position.x, vp.position.y, vp.position.z, vp.velocity.x, vp.velocity.y, vp.velocity.z); }
        println!();
    }
}
