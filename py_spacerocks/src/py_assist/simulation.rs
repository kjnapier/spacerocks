use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::assist::{PerturberCache, SpiceSimulation};
use std::sync::Arc;
   

use crate::PySpaceRock;
// use crate::py_spacerock::rockcollection::RockCollection;
use crate::py_time::time::PyTime;
use crate::py_assist::force::PyForce;
use crate::PySpiceKernel;


#[pyclass]
#[pyo3(name = "SpiceSimulation")]
pub struct PySpiceSimulation {
    pub inner: SpiceSimulation,
}

#[pymethods]
impl PySpiceSimulation {

    // #[new]
    // pub fn new() -> PyResult<Self> {
    //     match Simulation::new(&Time::now(), "J2000", "SSB") {
    //         Ok(sim) => Ok(PySimulation { inner: sim }),
    //         Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    //     }
    // }

    /// Create a new simulation with the giant planets.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, &'static str>` - The Simulation object.
    // #[classmethod]
    // pub fn giants(_cls: Py<PyType>, epoch: &PyTime, kernel: &PySpiceKernel) -> PyResult<Self> {
    //     match SpiceSimulation::giants(&epoch.inner, &kernel.inner) {
    //         Ok(sim) => Ok(PySpiceSimulation { inner: sim }),
    //         Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    //     }
    // }

    /// Create a new simulation with the JPL horizons perturbers.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, &'static str>` - The Simulation object.
    #[classmethod]
    pub fn horizons(_cls: Py<PyType>, epoch: &PyTime, kernel: &PySpiceKernel) -> PyResult<Self> {
        // need to Arc the kernel to pass it to the simulation
        match SpiceSimulation::horizons(&epoch.inner, &kernel.inner) {
            
            Ok(sim) => Ok(PySpiceSimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Add a SpaceRock to the simulation.
    ///
    /// # Arguments
    ///
    /// * `rock` - The SpaceRock to add.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn add(&mut self, rock: &PySpaceRock) -> PyResult<()> {
        match self.inner.add(rock.inner.clone()) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Integrate the simulation to a specific epoch.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate to.
    pub fn integrate(&mut self, py: Python<'_>, epoch: &PyTime, kernel: &PySpiceKernel) -> PyResult<()> {
        let inner = &mut self.inner;
        let k = &kernel.inner;
        let t = epoch.inner.clone();
        py.detach(|| inner.integrate(&t, k).map_err(|e| e.to_string()))
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }
    

    /// Step the simulation by one timestep.
    pub fn step(&mut self, kernel: &PySpiceKernel) -> PyResult<()> {
        self.inner.step(&kernel.inner).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Use a precomputed perturber ephemeris (a `PerturberCache`) for epochs it covers.
    /// Worth it for simulations with few particles that are integrated many times over the
    /// same span (e.g. orbit fitting); the cache can be shared between simulations.
    pub fn set_perturber_cache(&mut self, cache: PyRef<PyPerturberCache>) {
        self.inner.set_perturber_cache(cache.inner.clone());
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add.
    pub fn add_force(&mut self, force: PyRef<PyForce>) {
        self.inner.add_force(force.inner.clone());
    }

    /// Replace the force model, e.g. `sim.set_forces([Force.gr_eih(kernel, sources=11),
    /// Force.newtonian_gravity()])`. Forces are applied in order; put the smallest first.
    pub fn set_forces(&mut self, forces: Vec<PyRef<PyForce>>) {
        self.inner.set_forces(forces.iter().map(|f| f.inner.clone()).collect());
    }

    /// Use Newtonian point-mass gravity only.
    pub fn newtonian_only(&mut self) {
        self.inner.newtonian_only();
    }

    /// Add a variational particle for `parent` (a particle name): `dimension` is one of
    /// "x", "y", "z", "vx", "vy", "vz" (initial-state partials) or "A1", "A2", "A3"
    /// (non-gravitational parameter partials).
    pub fn add_variation(&mut self, dimension: &str, parent: &str) -> PyResult<()> {
        self.inner.add_variation(dimension, parent).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Add the six state variational particles of `parent`; with `nongrav=True` also the
    /// three for A1, A2, A3.
    #[pyo3(signature = (parent, nongrav = false))]
    pub fn add_full_variation(&mut self, parent: &str, nongrav: bool) -> PyResult<()> {
        let err = |e: Box<dyn std::error::Error>| pyo3::exceptions::PyValueError::new_err(e.to_string());
        self.inner.add_full_variation(parent).map_err(err)?;
        if nongrav {
            for a in ["A1", "A2", "A3"] {
                self.inner.add_variation(a, parent).map_err(err)?;
            }
        }
        Ok(())
    }

    /// States of the variational particles at the current epoch, shape (n, 6), in the order
    /// they were added. After `add_full_variation(name, nongrav=True)`, the transpose of the
    /// first six rows is the state transition matrix and the last three rows are the partials
    /// with respect to A1, A2, A3.
    #[getter]
    pub fn variational_particles<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        use numpy::IntoPyArray;
        let v = &self.inner.state.variational_particles;
        let flat: Vec<f64> = v
            .iter()
            .flat_map(|p| [p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z])
            .collect();
        numpy::ndarray::Array2::from_shape_vec((v.len(), 6), flat).unwrap().into_pyarray(py)
    }

    // /// Get a single particle from the simulation by name.
    // pub fn get_particle(&self, name: &str) -> PyResult<PySpaceRock> {
    //     let rock = self.inner.get_particle(name);
    //     match rock {
    //         Ok(r) => Ok(PySpaceRock { inner: r.clone() }),
    //         Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    //     }
    // }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    /// Length (days) of the last completed step.
    #[getter]
    pub fn last_timestep(&self) -> f64 {
        self.inner.integrator.last_timestep()
    }

    /// Epoch (TDB Julian date) at the end of the last completed step.
    #[getter]
    pub fn step_epoch(&self) -> f64 {
        let s = &self.inner.state;
        s.particles_1.first().map(|p| s.jd_ref + p.epoch).unwrap_or(s.epoch)
    }

    /// IAS15 tolerance (default 1e-9).
    #[getter]
    pub fn epsilon(&self) -> f64 {
        self.inner.integrator.epsilon()
    }

    #[setter]
    pub fn set_epsilon(&mut self, epsilon: f64) {
        self.inner.integrator.set_epsilon(epsilon);
    }

    /// Step-size criterion: "prs23" (default; Pham, Rein & Spiegel 2024), "global" (REBOUND's
    /// original criterion, which ASSIST uses), "per_particle" (the global criterion applied to
    /// each particle separately) or "individual".
    #[getter]
    pub fn adaptive_mode(&self) -> &'static str {
        self.inner.integrator.adaptive_mode().as_str()
    }

    #[setter]
    pub fn set_adaptive_mode(&mut self, mode: &str) -> PyResult<()> {
        let m = spacerocks::assist::AdaptiveMode::from_str(mode).map_err(pyo3::exceptions::PyValueError::new_err)?;
        self.inner.integrator.set_adaptive_mode(m);
        Ok(())
    }

    /// Round-off control in IAS15's sums: "kahan" (default; compensated positions and
    /// velocities), "full" (also the corrector's sums, as REBOUND/ASSIST do) or "none".
    #[getter]
    pub fn summation(&self) -> &'static str {
        self.inner.integrator.summation().as_str()
    }

    #[setter]
    pub fn set_summation(&mut self, mode: &str) -> PyResult<()> {
        let m = spacerocks::assist::Summation::from_str(mode).map_err(pyo3::exceptions::PyValueError::new_err)?;
        self.inner.integrator.set_summation(m);
        Ok(())
    }

    /// Smallest allowed step (days; default 1e-8).
    #[getter]
    pub fn min_timestep(&self) -> f64 {
        self.inner.integrator.min_timestep()
    }

    #[setter]
    pub fn set_min_timestep(&mut self, min_timestep: f64) {
        self.inner.integrator.set_min_timestep(min_timestep);
    }

    #[getter]
    pub fn particles(&self) -> Vec<PySpaceRock> {
        self.inner.state.particles.iter().map(|rock| PySpaceRock { inner: rock.clone() }).collect()
    }

}

/// Chebyshev fits of the `SpiceSimulation.horizons` perturbers over a time span.
///
/// `PerturberCache(kernel, start, end)` fits the barycentric states of the Sun, planets, Moon,
/// Pluto and 16 asteroids between two epochs (to ~1e-14 AU); `SpiceSimulation.set_perturber_cache`
/// then reads perturber states from the fits instead of the kernel, about twice as fast.
#[pyclass]
#[pyo3(name = "PerturberCache")]
pub struct PyPerturberCache {
    pub inner: Arc<PerturberCache>,
}

#[pymethods]
impl PyPerturberCache {
    #[new]
    fn new(py: Python<'_>, kernel: &PySpiceKernel, start: &PyTime, end: &PyTime) -> PyResult<Self> {
        let ids = SpiceSimulation::horizons_body_ids().map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let (a, b) = (start.inner.tdb().jd(), end.inner.tdb().jd());
        let k = &kernel.inner;
        let cache = py
            .detach(|| PerturberCache::build(k, &ids, a.min(b), a.max(b)).map_err(|e| e.to_string()))
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        Ok(PyPerturberCache { inner: Arc::new(cache) })
    }

    /// Cached span as (start, end) TDB Julian dates.
    #[getter]
    fn span(&self) -> (f64, f64) {
        self.inner.span()
    }

    fn __repr__(&self) -> String {
        let (a, b) = self.inner.span();
        format!("PerturberCache({} bodies, JD {} to {} TDB)", self.inner.ids().len(), a, b)
    }
}
