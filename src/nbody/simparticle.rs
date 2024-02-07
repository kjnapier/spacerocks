use crate::spacerock::SpaceRock;
use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct SimParticle {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
    pub epoch: f64,
    pub mass: f64,
    pub name: String,
}

impl SimParticle {

    pub fn from_spacerock(spacerock: &SpaceRock) -> Self {

        let m;
        if let Some(mass) = spacerock.mass {
            m = mass;
        } else {
            m = 0.0;
        }
           
        SimParticle {
            name: spacerock.name.clone(),
            position: spacerock.position,
            velocity: spacerock.velocity,
            acceleration: Vector3::new(0.0, 0.0, 0.0),
            epoch: spacerock.epoch.epoch,
            mass: m,
        }
    }

}