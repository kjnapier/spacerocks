use spacerocks::{SpaceRock, Time, SpiceKernel, Observatory, Simulation};
use spacerocks::SpiceSimulation;

use std::sync::Arc;

use nalgebra::Vector3;
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let spice_root = "/Users/kjnapier/data/spice";
    // let spice_root = "/home/kevin/data/spice/";

    // load spice kernels
    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", spice_root).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", spice_root).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", spice_root).as_str())?;

    let mut epoch = Time::new(2460762.549988426, "TDB", "JD")?;
    let mut rock = SpaceRock::from_xyz("holman", 2.963305899720348, -1.627586306680811, -0.7799786968810375, 
                                        0.004951381894813546, 0.006677060249157604, 0.002540378471598749, epoch.clone(), "J2000", "SSB")?;

    // Need to wrap the kernel in an Arc, since we don't want to clone it or give ownership to the simulation.
    // The simulation gets a pointer to the kernel, and the kernel, and the kernel can be shared between multiple simulations.
    // Since the simulations are not modifying the kernel, they can share it safely.
    let mut sim = SpiceSimulation::horizons(&epoch, &kernel)?;
    sim.add(rock)?;


    let mut times = Vec::new();
    for idx in 0..1000 {
        let t = epoch.clone() + (idx as f64);
        times.push(t);
    }

    let start = std::time::Instant::now();
    for t in &times {
        sim.integrate(t, &kernel)?;
    }
    let elapsed = start.elapsed();
    println!("Time to integrate: {:?}", elapsed);
    
  
    Ok(())
}