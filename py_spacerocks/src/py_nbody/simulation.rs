use std::sync::{Arc, Mutex};

use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::nbody::Simulation;
use spacerocks::coordinates::{ReferencePlane, Origin};
use spacerocks::Time;
use numpy::{IntoPyArray, PyArray1, PyArray2};   

use crate::PySpaceRock;
// use crate::py_spacerock::rockcollection::RockCollection;
use crate::py_time::time::PyTime;
use crate::py_nbody::integrator::PyIntegrator;
use crate::py_nbody::force::PyForce;
use crate::py_coordinates::origin::PyOrigin;
use crate::PySpiceKernel;

#[pyclass]
#[pyo3(name = "Simulation")]
pub struct PySimulation {
    pub inner: Simulation,
}

#[pymethods]
impl PySimulation {

    #[new]
    pub fn new() -> PyResult<Self> {
        match Simulation::new(&Time::now(), "J2000", "SSB") {
            Ok(sim) => Ok(PySimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

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
    #[classmethod]
    pub fn giants(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::giants(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Create a new simulation with the planets.
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
    pub fn planets(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::planets(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

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
    pub fn horizons(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::horizons(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
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

    /// Remove a SpaceRock from the simulation.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the SpaceRock to remove.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn remove(&mut self, name: &str) -> PyResult<()> {
        match self.inner.remove(name) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Integrate the simulation to a specific epoch.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate to.
    pub fn integrate(&mut self, py: Python<'_>, epoch: &PyTime) -> PyResult<()> {
        let epoch = epoch.inner.clone();
        run_interruptible(py, &mut self.inner, move |sim| sim.integrate(&epoch))
    }

    /// Integrate with no end, until the monitor page's Stop (or an end set there), or Ctrl-C
    /// (which raises KeyboardInterrupt). Goes backwards with a negative timestep.
    pub fn run(&mut self, py: Python<'_>) -> PyResult<()> {
        run_interruptible(py, &mut self.inner, |sim| sim.run())
    }

    /// Whether solar general relativity is among the forces.
    #[getter]
    pub fn solar_gr(&self) -> bool {
        self.inner.solar_gr()
    }

    /// Add solar general relativity to the forces (True), or take it out (False).
    pub fn set_solar_gr(&mut self, on: bool) {
        self.inner.set_solar_gr(on);
    }

    /// Turn the gravity of the massive particle `name` off (it then moves as a test particle)
    /// or back on. The central body can't be turned off.
    pub fn set_perturber(&mut self, name: &str, on: bool) -> PyResult<()> {
        self.inner.set_perturber(name, on).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }

    /// Every massive particle but the central body, as (name, gravity on).
    #[getter]
    pub fn perturbers(&self) -> Vec<(String, bool)> {
        self.inner.perturbers()
    }

    /// Integrate to an epoch with the integrator's own steps, then interpolate to it.
    ///
    /// Much faster than `integrate` for many output epochs. Falls back to `integrate` for
    /// integrators without dense output (only IAS15 has it).
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate or interpolate to.
    pub fn integrate_or_interpolate(&mut self, py: Python<'_>, epoch: &PyTime) -> PyResult<()> {
        let epoch = epoch.inner.clone();
        run_interruptible(py, &mut self.inner, move |sim| sim.integrate_or_interpolate(&epoch))
    }
    
    /// Serve a live progress page for this simulation and return its URL.
    ///
    /// The page shows how far `integrate`/`integrate_or_interpolate` has got towards its target
    /// epoch, the current epoch and timestep, steps per second, the time left, the energy error
    /// of the massive bodies and a top-down view of the particles. It listens on 127.0.0.1 only.
    ///
    /// # Arguments
    ///
    /// * `port` - The port to serve on; 0 (the default) picks a free one.
    /// * `open_browser` - Open the page in the default web browser.
    /// * `paused` - Start paused, so the next integration waits for Resume or Step on the page.
    /// * `track` - Name of a particle whose orbital elements the page follows (clicking a
    ///   particle on the page does the same).
    ///
    /// The page's Pause holds the simulation before its next step, Step lets one step through,
    /// and the integrator and timestep can be switched there. While paused, the call that is
    /// integrating waits until the page resumes it, or until Ctrl-C.
    #[pyo3(signature = (port=0, open_browser=true, paused=false, track=None))]
    pub fn monitor(&mut self, py: Python<'_>, port: u16, open_browser: bool, paused: bool, track: Option<&str>) -> PyResult<String> {
        let url = self.inner.monitor(port).map_err(|e| PyErr::new::<pyo3::exceptions::PyOSError, _>(e.to_string()))?;
        if let Some(monitor) = self.inner.progress_monitor() {
            if paused {
                monitor.set_paused(true);
            }
            if track.is_some() {
                monitor.track(track);
            }
        }
        if open_browser {
            py.import("webbrowser")?.call_method1("open", (url.as_str(),))?;
        }
        Ok(url)
    }

    /// Measure the progress page's progress from the current epoch up to `epoch`, across however
    /// many calls it takes (e.g. a loop of `integrate_or_interpolate` over output epochs),
    /// instead of over each call.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch the whole run ends at.
    pub fn set_monitor_target(&mut self, epoch: &PyTime) {
        self.inner.set_monitor_target(&epoch.inner);
    }

    /// Detach the progress page and stop its server.
    pub fn stop_monitor(&mut self) {
        self.inner.stop_monitor();
    }

    /// Follow the orbital elements of the particle called `name` on the progress page (None stops).
    #[pyo3(signature = (name=None))]
    pub fn monitor_track(&self, name: Option<&str>) {
        if let Some(monitor) = self.inner.progress_monitor() {
            monitor.track(name);
        }
    }

    /// Pause or resume the simulation, as the progress page's Pause button does.
    pub fn set_monitor_paused(&self, paused: bool) {
        if let Some(monitor) = self.inner.progress_monitor() {
            monitor.set_paused(paused);
        }
    }

    /// The progress page's URL, or None.
    #[getter]
    pub fn monitor_url(&self) -> Option<String> {
        self.inner.monitor_url().map(str::to_string)
    }

    /// Step the simulation by one timestep.
    pub fn step(&mut self, py: Python<'_>) -> PyResult<()> {
        run_interruptible(py, &mut self.inner, |sim| sim.step())
    }

    /// Step the simulation by `n` timesteps: the same result as calling `step` `n` times, but
    /// WisdomHolman and Trace then split test particles between threads.
    pub fn steps(&mut self, py: Python<'_>, n: usize) -> PyResult<()> {
        run_interruptible(py, &mut self.inner, move |sim| sim.steps(n))
    }

    /// Move the simulation to the center of mass.
    pub fn move_to_center_of_mass(&mut self) -> PyResult<()> {
        match self.inner.move_to_center_of_mass() {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Change the origin of the simulation. The new origin must be present in the simulation.
    ///
    /// # Arguments
    ///
    /// * `origin` - The new origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn change_origin(&mut self, origin: &str) -> PyResult<()> {
        match self.inner.change_origin(origin) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Change the reference plane of the simulation.
    ///
    /// # Arguments
    ///
    /// * `reference_plane` - The new reference plane of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn set_epoch(&mut self, epoch: &PyTime) {
        self.inner.set_epoch(&epoch.inner);
    }

    /// Set the reference plane of the simulation.
    ///
    /// # Arguments
    ///
    /// * `reference_plane` - The new reference plane of the simulation.
    pub fn set_reference_plane(&mut self, reference_plane: &str) {
        let reference_plane = ReferencePlane::from_str(reference_plane).unwrap();
        self.inner.reference_plane = reference_plane;
    }

    /// Set the origin of the simulation.
    ///
    /// # Arguments
    ///
    /// * `origin` - The new origin of the simulation.
    pub fn set_origin(&mut self, origin: &str) -> PyResult<()> {
        match Origin::from_str(origin) {
            Ok(o) => {
                self.inner.origin = o;
                Ok(())
            },
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Set the integrator.
    ///
    /// # Arguments
    ///
    /// * `integrator` - The new integrator.
    pub fn set_integrator(&mut self, integrator: PyRef<PyIntegrator>) {
        self.inner.integrator = integrator.inner.clone();
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add.
    pub fn add_force(&mut self, force: PyRef<PyForce>) {
        self.inner.add_force(force.inner.clone());
    }

    /// Calculate the total energy of the simulation.
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }

    /// States of the particles as an `(n, 6)` array of `x, y, z` (AU) and `vx, vy, vz`
    /// (AU/day), massive particles first. A copy.
    #[getter]
    pub fn states<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let n = self.inner.len();
        let flat = self.inner.states().as_flattened().to_vec();
        numpy::ndarray::Array2::from_shape_vec((n, 6), flat).unwrap().into_pyarray(py)
    }

    /// Masses of the particles (solar masses), in the order of `states`.
    #[getter]
    pub fn masses<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.masses().to_vec().into_pyarray(py)
    }

    /// Names of the particles, in the order of `states`.
    #[getter]
    pub fn names(&self) -> Vec<String> {
        self.inner.names().to_vec()
    }

    /// Every particle as a SpaceRock, in the order of `states`.
    pub fn particles(&self) -> Vec<PySpaceRock> {
        self.inner.particles().into_iter().map(|r| PySpaceRock { inner: r }).collect()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Get a single particle from the simulation by name.
    pub fn get_particle(&self, name: &str) -> PyResult<PySpaceRock> {
        let rock = self.inner.get_particle(name);
        match rock {
            Ok(r) => Ok(PySpaceRock { inner: r.clone() }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
       
    }


    // #[getter]
    // pub fn particles(&self) -> RockCollection {
    //     RockCollection { rocks: self.inner.particles.clone(), name_hash_map: self.inner.particle_index_map.clone() }
    // }


    #[getter]
    pub fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch() }
    }

    #[getter]
    pub fn reference_plane(&self) -> String {
        self.inner.reference_plane.to_string()
    }


    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    #[getter]
    pub fn origin(&self) -> PyOrigin {
        let o = self.inner.origin.clone();
        PyOrigin { inner: o }
    }


    // repr
    pub fn __repr__(&self) -> String {
        // format!("Simulation:\n");
        // format!("Epoch: {}\n", self.inner.epoch);
        // format!("Reference Plane: {}\n", self.inner.reference_plane);
        // format!("Origin: {}\n", self.inner.origin);
        // format!("Timestep: {}\n", self.inner.integrator.timestep());
        // format!("Forces: {:?}\n", self.inner.forces);
        // format!("Integrator: {:?}\n", self.inner.integrator);
        // format!("Particles: {}\n", self.inner.particles.len());

        let mut s = String::new();
        s.push_str("Simulation:\n");
        s.push_str(&format!("    Epoch: {}\n", self.inner.epoch()));
        s.push_str(&format!("    Reference Plane: {}\n", self.inner.reference_plane));
        s.push_str(&format!("    Origin: {}\n", self.inner.origin));
        s.push_str(&format!("    Timestep: {}\n", self.inner.integrator.timestep()));
        s
    }

    #[getter]
    pub fn integrator(&self) -> PyIntegrator {
        PyIntegrator { inner: self.inner.integrator.clone() }
    }

}

/// Run `f` on the simulation without holding the GIL, so other Python threads keep running, and
/// stop it between steps when Python has a signal to handle (Ctrl-C raises KeyboardInterrupt,
/// with the simulation left at the last step it finished).
fn run_interruptible<F>(py: Python<'_>, sim: &mut Simulation, f: F) -> PyResult<()>
where
    F: FnOnce(&mut Simulation) + Send,
{
    let error: Arc<Mutex<Option<PyErr>>> = Arc::new(Mutex::new(None));
    let seen = error.clone();
    sim.set_interrupt(Some(Arc::new(move || {
        Python::attach(|py| match py.check_signals() {
            Ok(()) => false,
            Err(e) => {
                *seen.lock().unwrap() = Some(e);
                true
            }
        })
    })));
    py.detach(|| f(&mut *sim));
    sim.set_interrupt(None);
    let caught = error.lock().unwrap().take();
    match caught {
        Some(e) => Err(e),
        None => Ok(()),
    }
}
