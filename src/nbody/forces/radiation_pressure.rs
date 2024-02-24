use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;


pub struct RadiationPressure;

impl Force for RadiationPressure {

    fn apply(&self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>) {
        _;
    }

}