use spacerocks::{SpaceRock, Time, SpiceKernel, Observatory, Simulation};
use spacerocks::SpiceSimulation;

use std::sync::Arc;

use nalgebra::Vector3;

use plotly::{Plot, Scatter};

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


    let mut perturbed_rock = SpaceRock::from_xyz("holman_dx", 2.963305899720348 + 1.0e-8, -1.627586306680811, -0.7799786968810375, 
                                        0.004951381894813546, 0.006677060249157604, 0.002540378471598749, epoch.clone(), "J2000", "SSB")?;

    // Need to wrap the kernel in an Arc, since we don't want to clone it or give ownership to the simulation.
    // The simulation gets a pointer to the kernel, and the kernel, and the kernel can be shared between multiple simulations.
    // Since the simulations are not modifying the kernel, they can share it safely.
    let mut sim = SpiceSimulation::horizons(&epoch, &kernel)?;
    sim.add(rock)?;
    sim.add(perturbed_rock)?;

    sim.add_variation("x", "holman");


    // sim.step(&kernel);
    
    let mut times = Vec::new();
    for idx in 0..20000 {
        let t = epoch.clone() + (idx as f64);
        times.push(t);
    }

    let mut shadow_dx = Vec::new();
    let mut analytic_dx = Vec::new();

    let start = std::time::Instant::now();
    for t in &times {
        sim.integrate(t, &kernel)?;
        let r1 = sim.state.particles[1].position.x - sim.state.particles[0].position.x;
        shadow_dx.push(r1);

        let r2 = sim.state.variational_particles[0].position.x * 1.0e-8;
        analytic_dx.push(r2);
        println!("{:?}", sim.state.variational_particles[0].position);
    }
    let elapsed = start.elapsed();
    println!("Time to integrate: {:?}", elapsed);

    let au_km = 149597870.700;
    // Convert to km
    shadow_dx = shadow_dx.iter().map(|x| x * au_km).collect::<Vec<f64>>();
    analytic_dx = analytic_dx.iter().map(|x| x * au_km).collect::<Vec<f64>>();

    let diffs = shadow_dx.iter().zip(analytic_dx.iter()).map(|(x, y)| x - y).collect::<Vec<f64>>();
    // multiply by 1e6 to get mm
    let diffs = diffs.iter().map(|x| x * 1.0e6).collect::<Vec<f64>>();


    let ts = times.iter().map(|t| t.tdb().jd()).collect::<Vec<f64>>();

    // Create a scatter trace
    let trace = Scatter::new(ts.clone(), shadow_dx)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Shadow Particle");

    // Create a scatter trace for the analytic solution
    let trace_analytic = Scatter::new(ts.clone(), analytic_dx)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Analytic Solution");

    // Create and render the plot
    let mut plot = Plot::new();
    plot.add_trace(trace);

    // open the plot in a new window
    plot.add_trace(trace_analytic);
    

    plot.show();
    

    let mut plot = Plot::new();
    let trace = Scatter::new(ts.clone(), diffs)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Shadow Particle");
    plot.add_trace(trace);
    plot.show();



  
    Ok(())
}