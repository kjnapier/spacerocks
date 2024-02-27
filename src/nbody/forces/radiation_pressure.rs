use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;

#[derive(Debug, Clone, Copy)]
pub struct RadiationPressure;

impl Force for RadiationPressure {

    fn apply(&self, entities: &mut Vec<SpaceRock>) {
        _;
    }

}