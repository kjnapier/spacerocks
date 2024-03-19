use pyo3::prelude::*;

use spacerocks::orbfit::gauss;
use spacerocks::observing::detection::Detection;

use crate::py_spacerock::spacerock::PySpaceRock;
use crate::py_observing::detection::PyDetection;


// #[pyfunction]
// fn gauss_fit(dets:  

    // takes a list of three pydetections 

#[pyfunction]
pub fn gauss_fit(dets: Vec<PyRef<PyDetection>>, min_distance: f64) -> PyResult<Vec<PySpaceRock>> {
    
    let mut detections: Vec<Detection> = Vec::new();
    for det in dets {
        detections.push(det.inner.clone());
    }

    let mut dets: Vec<&Detection> = Vec::new();
    for det in detections.iter() {
        dets.push(det);
    }

    let result = gauss::gauss_fit(&dets, min_distance);

    match result {
        Some(rocks) => {
            let mut pyrocks = Vec::new();
            for rock in rocks {
                pyrocks.push(PySpaceRock { inner: rock });
            }
            Ok(pyrocks)
        },
        None => Ok(Vec::new())
    }
}