use spacerocks::spacerock::SpaceRock;
use spacerocks::properties::Properties;
use spacerocks::observatory::Observatory;
use spacerocks::nbody::Simulation;
use spacerocks::time::Time;

use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let root = "/Users/kjnapier/Desktop/research/data/spice";
    let outdir = "/Users/kjnapier/Desktop";

    spice::furnsh(&format!("{root}/latest_leapseconds.tls", root=root));
    spice::furnsh(&format!("{root}/de440s.bsp", root=root));

    // pull the rock from JPL Horizons

    let mut arrokoth = SpaceRock::from_horizons("Arrokoth")?;
    
    let properties = Properties { H: Some(11.06), Gslope: Some(0.15), ..Default::default() };
    arrokoth.properties = Some(properties);

    let t0 = arrokoth.epoch.clone();

    // pull the planets from spice
    let sun = SpaceRock::from_spice("Sun", &t0);
    let jupiter = SpaceRock::from_spice("Jupiter Barycenter", &t0);
    let saturn = SpaceRock::from_spice("Saturn Barycenter", &t0);
    let uranus = SpaceRock::from_spice("Uranus Barycenter", &t0);
    let neptune = SpaceRock::from_spice("Neptune Barycenter", &t0);

    let mut sim = Simulation::new();
    sim.time = t0.epoch;
    sim.timestep = 0.001;

    let perturbers = vec![sun, jupiter, saturn, uranus, neptune];
    for mut perturber in perturbers {
        perturber.change_frame("ECLIPJ2000");
        sim.add(perturber);
    }

    arrokoth.change_frame("ECLIPJ2000");
    sim.add(arrokoth.clone());
    sim.move_to_center_of_mass();

    let w84 = Observatory::from_coordinates(-30.00293494202556, -70.80642, 2207.0);

    let mut epochs: Vec<Time> = Vec::new();
    for idx in 0..50_000 {
        let t = Time::new(t0.epoch + idx as f64 * 0.1, "tdb", "jd");
        epochs.push(t);
    }

    let mut file = File::create(&format!("{outdir}/positions.csv", outdir=outdir)).unwrap();
    file.write_all(b"objid,epoch,ra,dec,mag\n").unwrap();

    for epoch in epochs {
        let observer = w84.at(&epoch);
        let sun = SpaceRock::from_spice("Sun", &epoch);
        sim.integrate(epoch.epoch);
        let test_particle = sim.get_particle("Arrokoth").ok_or("Arrokoth not found")?;
        let obs = test_particle.clone().observe(&observer);
        let mag = obs.mag(&sun);
        file.write_all(format!("{objid}, {epoch}, {ra}, {dec}, {mag}\n", 
                                objid=test_particle.name, 
                                epoch=test_particle.epoch.epoch, 
                                ra=obs.ra, 
                                dec=obs.dec,
                                mag=mag
                                ).as_bytes())?;
    }

    arrokoth.analytic_propagate(&t0);
    return Ok(());

}