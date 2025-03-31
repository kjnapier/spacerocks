use spacerocks::{SpaceRock, Time, SpiceKernel, Observatory, Simulation};
use spacerocks::SpiceSimulation;

use std::sync::Arc;

use nalgebra::Vector3;
fn main() -> Result<(), Box<dyn std::error::Error>> {

    // let spice_root = "/Users/kjnapier/data/spice";
    let spice_root = "/home/kevin/data/spice/";

    // load spice kernels
    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", spice_root).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", spice_root).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", spice_root).as_str())?;

    let epoch = Time::new(2451545.0, "tdb", "jd")?;
    let mut rock = SpaceRock::from_horizons("holman", &epoch, "J2000", "SSB")?;

    // let mut epoch = Time::new(2460762.549988426, "UTC", "JD")?;
    // let mut rock = SpaceRock::from_xyz("holman", 2.963305899720348, -1.627586306680811, -0.7799786968810375, 
    //                                     0.004951381894813546, 0.006677060249157604, 0.002540378471598749, epoch.clone(), "J2000", "SSB")?;

    // Need to wrap the kernel in an Arc, since we don't want to clone it or give ownership to the simulation.
    // The simulation gets a pointer to the kernel, and the kernel, and the kernel can be shared between multiple simulations.
    // Since the simulations are not modifying the kernel, they can share it safely.
    let kernel = Arc::new(kernel);
    let mut sim = SpiceSimulation::horizons(&epoch, "J2000", "SSB", &kernel)?;


    let start = std::time::Instant::now();
    for idx in 0..10000 {
        kernel.get_barycentric_states(&sim.state.spice_bodies, sim.state.epoch + (idx as f64) * 1.0)?;
    }
    let elapsed = start.elapsed();
    let time_per_step = elapsed.as_secs_f64() / 10000.0;
    println!("Time per step: {:?} us", time_per_step * 1_000_000.0);
    // let mut sim = SpiceSimulation::horizons(&epoch, "J2000", "SSB", &kernel)?;

    sim.add(rock)?;

    let start = std::time::Instant::now();
    for _ in 0..1000 {
        sim.step();
    }
    let elapsed = start.elapsed();
    let time_per_step = elapsed.as_secs_f64() / 1000.0;
    println!("Time per step: {:?} us", time_per_step * 1_000_000.0);


    // sim.step();
    // sim.step();
    // sim.step();
    // sim.step();
    // println!("sim: {:?}", sim.state.particles_1[0].position);

    
    // let mut rock = SpaceRock::from_horizons("holman", &Time::new(sim.state.epoch, "tdb", "jd")?, "J2000", "SSB")?;
    // println!("jpl: {:?}", rock.position);


    // let t2 = epoch.clone() + 1000.0;
    // sim.integrate(&t2);
    // println!("sim: {:?}", sim.state.particles[0].position);


    // let mut rock = SpaceRock::from_horizons("holman", &t2, "J2000", "SSB")?;
    // println!("jpl: {:?}", rock.position);

    
    // for _ in 0..1000 {
    //     sim.step();
    // }
    
    // println!("Simulation created successfully");
    // println!("sim: {:?}", sim.state.particles[0].position);


    // let mut rock = SpaceRock::from_horizons("holman", &Time::new(sim.state.epoch, "tdb", "jd")?, "J2000", "SSB")?;
    // println!("jpl: {:?}", rock.position);


    

    // let start = std::time::Instant::now();
    // kernel.get_barycentric_states(&sim.state.spice_bodies, sim.state.epoch)?;
    // let elapsed = start.elapsed();
    // println!("Time to get barycentric states: {:?}", elapsed);

   
    


    // println!("Simulation created successfully");
    // println!("sim: {:?}", sim.state.spice_particles);

    // // observer not working rn, so commenting out for complilation
    // let f51 = Observatory::from_obscode("F51")?;

    // let observer = f51.at(&epoch, "J2000", "SSB")?;
    // println!("{:?}", observer);

    // let mut arrokoth = SpaceRock::from_horizons("Arrokoth", &epoch, "J2000", "SSB")?;
    // println!("{}", arrokoth);

    // let mut sim = Simulation::horizons(&epoch, "J2000", "SSB")?;
    // println!("Sim is working");
    // sim.add(arrokoth.clone())?;

    // println!("{}", sim.epoch);

    // let future = epoch + 365.25 * 10.0;

    // println!("Arrkoth position and velocity at epoch: {}, {}", &arrokoth.position, &arrokoth.velocity);

    // arrokoth.analytic_propagate(&future)?;

    // println!("Arrkoth position and velocity at future: {}, {}", &arrokoth.position, &arrokoth.velocity);

    // let observer = f51.at(&epoch)?;

    // let mut arrokoth = SpaceRock::from_horizons("Arrokoth", &epoch, "j2000", "ssb")?;
    // let observation = arrokoth.observe(&observer)?;
    // let mut arrokoth = SpaceRock::from_horizons("Arrokoth", &epoch, "j2000", "ssb")?;
    // let observation = arrokoth.observe(&observer)?;

    // println!("{}", observation);
    // println!("{}", observation);
    


    // // load arrokoth from spice
    // let rock = SpaceRock::from_horizons("Arrokoth", &epoch, "j2000", "ssb")?;
    // // println!("{}", rock);

    // create a simulation and add arrokoth
    // let mut sim = Simulation::giants(&epoch, "ECLIPJ2000", "SSB")?;
    // sim.add(rock)?; 
    // sim.move_to_center_of_mass()?;

    // let dt = 10.0;
    // let t_total = 365.25 * 10_000.0;
    // let n_epochs = (t_total / dt) as usize;

    // let mut positions: HashMap<String, Vec<Vector3<f64>>> = HashMap::new();
    // for i in 0..n_epochs {
    //     let t = epoch.clone() + (i as f64) * dt;
    //     sim.integrate(&t);
    //     for rock in &sim.particles {
    //         let pos = rock.position;
    //         let name = rock.name.clone();
    //         if positions.contains_key(&name) {
    //             positions.get_mut(&name).unwrap().push(pos);
    //         } else {
    //             positions.insert(name.to_string(), vec![pos]);
    //         }
    //     }
    // }
    
    // let caption = format!("{} year simulation", (t_total / 365.25) as usize);
    // let root = BitMapBackend::new("plot.png", (800, 800)).into_drawing_area();
    // root.fill(&WHITE)?;
    // let mut chart = ChartBuilder::on(&root)
    //     .caption(caption, ("sans-serif", 30).into_font())
    //     .margin(5)
    //     .x_label_area_size(40)
    //     .y_label_area_size(40)
    //     .build_ranged(-50.0..50.0, -50.0..50.0)?;
    // chart.configure_mesh().draw()?;

    // for (name, pos) in positions.iter() {
    //     let x: Vec<f64> = pos.iter().map(|v| v.x).collect();
    //     let y: Vec<f64> = pos.iter().map(|v| v.y).collect();
    //     chart.draw_series(
    //         x.iter().zip(y.iter()).map(|(x, y)| {
    //             Circle::new((*x, *y), 1, &BLACK)
    //         })
    //     )?;
    // }
   
    Ok(())
}