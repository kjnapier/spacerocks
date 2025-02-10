use crate::SpaceRock;

use nalgebra::{DVector, DMatrix};

pub struct AnalyticOrbitFitter {}

impl Model for AnalyticOrbitFitter {
    type Jacobian = Vec<Vec<f64>>; // In practice, you might use a dedicated matrix type.

    fn residuals(&self, rock: &mut SpaceRock, detections: &Vec<&Observation>) -> Vec<f64> {
        let mut res = Vec::new();
        for detection in detections.iter() {

            let observed_parameters = match detection.observation_type {
                ObservationType::Astrometry { ra, dec } => DVector::from_vec(vec![ra, dec]),
                ObservationType::Streak { ra, dec, ra_rate, dec_rate } => DVector::from_vec(vec![ra, dec, ra_rate, dec_rate]),
                _ => return Err("Observation type not supported".into()),
            };
    
            // somehow get the rock to the epoch of the detection
            rock.analytic_propagate(&detection.epoch);
    
            // calculate the model observations
            let astro = rock.observe(&detection.observer)?;
            let model_parameters = match detection.observation_type {
                ObservationType::Astrometry { ra, dec } => DVector::from_vec(vec![astro.ra(), astro.dec()]),
                ObservationType::Streak { ra, dec, ra_rate, dec_rate } => DVector::from_vec(vec![astro.ra(), astro.dec(), astro.ra_rate().expect("Must have ra rate"), astro.dec_rate().expect("Must have ra rate")]),
                _ => return Err("Observation type not supported".into()),
            };
    
            let d = observed_parameters - model_parameters;
            let m = &d.transpose() * &detection.inverse_covariance.clone().unwrap() * &d;
            residuals[idx] = m[0];
        }
        Ok(res)
    }

    // fn jacobian(&self, params: &Vec<f64>, detections: &Vec<&Observation>) -> Self::Jacobian {
    // }

    // fn residuals_and_derivatives(&self, params: &Vec<f64>, detections: &Vec<&Observation>) -> (Vec<f64>, Self::Jacobian) {
    // }

    // fn chisq(&self, params: &Vec<f64>, detections: &Vec<&Observation>) -> f64 {
    // }
}



