
use crate::SpaceRock;
use crate::Simulation;

#[derive(PartialEq, Debug, Clone)]
pub struct Freezer {
    pub timestep: f64,
}

impl Freezer {
    pub fn new(timestep: f64) -> Freezer {
        Freezer { timestep }
    }

    pub fn step(&self, sim: &mut Simulation) {
        for particle in &mut *sim.particles {
            particle.epoch += self.timestep;
        }

        for mut perturber in &mut *sim.perturbers {
            perturber.epoch += self.timestep;
        }
    }

}