// use pyo3::prelude::*;


// #[pyclass]
// #[pyo3(name = "Integrator")]
// pub struct PyIntegrator {
//     pub inner: crate::nbody::integrators::Integrator,
// }

// #[pymethods]
// impl PyIntegrator {

//     #[new]
//     pub fn new() -> Self {
//         PyIntegrator { inner: crate::nbody::integrators::Integrator::Leapfrog(crate::nbody::leapfrog::Leapfrog::new(20.0)) }
//     }

//     pub fn set_timestep(&mut self, timestep: f64) {
//         self.inner.set_timestep(timestep);
//     }

//     pub fn step(&self, perturbers: &mut Vec<crate::spacerock::SpaceRock>, particles: &mut Vec<crate::spacerock::SpaceRock>) {
//         self.inner.step(perturbers, particles);
//     }

//     pub fn threaded_step(&self, perturbers: &mut Vec<crate::spacerock::SpaceRock>, particles: &mut Vec<crate::spacerock::SpaceRock>) {
//         self.inner.threaded_step(perturbers, particles);
//     }

//     pub fn timestep(&self) -> f64 {
//         self.inner.timestep()
//     }

//     pub fn set_timestep(&mut self, timestep: f64) {
//         self.inner.set_timestep(timestep);
//     }
// }