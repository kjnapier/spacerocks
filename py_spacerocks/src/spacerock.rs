use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::SpaceRock;

use nalgebra::Vector3;

use crate::py_time::time::PyTime;
use crate::py_coordinates::origin::PyOrigin;
use crate::py_observing::observer::{PyObserver};
use crate::py_observing::observation::{PyObservation};
use crate::PySpiceKernel;

#[pyclass]
#[pyo3(name = "SpaceRock")]
pub struct PySpaceRock {
    pub inner: SpaceRock,
}

#[pymethods]
impl PySpaceRock {

    #[classmethod]
    #[pyo3(signature = (name, epoch, reference_plane="ECLIPJ2000", origin="SSB"))]
    fn from_horizons(_cls: Py<PyType>, name: &str, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str) -> PyResult<Self> {
        let rock = SpaceRock::from_horizons(name, &epoch.inner, reference_plane, origin);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Horizons for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
   //  #[pyo3(signature = (name, epoch, reference_plane="ECLIPJ2000", origin="SSB"))]
    fn from_spice(_cls: Py<PyType>, name: &str, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        let rock = SpaceRock::from_spice(name, &epoch.inner, reference_plane, origin, &kernel.inner);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Spice for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    #[pyo3(signature = (epoch, reference_plane="ECLIPJ2000", origin="SSB"))]
    fn random(_cls: Py<PyType>, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str) -> PyResult<Self> {
        let rock = SpaceRock::random(&epoch.inner, reference_plane, origin);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Failed to create random SpaceRock"));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    fn from_spherical(_cls: Py<PyType>, name: &str, phi: f64, theta: f64, r: f64, vr: f64, vo: f64, psi: f64, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str) -> PyResult<Self> {
        let rock = SpaceRock::from_spherical(name, phi, theta, r, vr, vo, psi, epoch.inner.clone(), reference_plane, origin);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Spherical for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    fn from_xyz(_cls: Py<PyType>, name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str) -> PyResult<Self> {
        let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, epoch.inner.clone(), reference_plane, origin);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from XYZ for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    fn from_kepler(_cls: Py<PyType>, name: &str, q: f64, e: f64, inc: f64, arg: f64, node: f64, true_anomaly: f64, epoch: PyRef<PyTime>, reference_plane: &str, origin: &str) -> PyResult<Self> {
        let rock = SpaceRock::from_kepler(name, q, e, inc, arg, node, true_anomaly, epoch.inner.clone(), reference_plane, origin);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Keplerian for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }


    fn observe(&mut self, observer: &PyObserver) -> PyResult<PyObservation> {
        // if observer.inner.frame != ReferencePlane::J2000 {
        //     return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks.")));
        // }

        match self.inner.observe(&observer.inner) {
            Ok(obs) => Ok(PyObservation { inner: obs }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to observe rock: {}", e))),
        }
    }

    fn propagate(&mut self, epoch: PyRef<PyTime>, kernel: &PySpiceKernel) -> PyResult<()> {
        let ker = &kernel.inner;
        match self.inner.propagate(&epoch.inner, &ker) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to propagate rock: {}", e))),
        }
    }

    fn analytic_propagate(&mut self, epoch: PyRef<PyTime>) -> PyResult<()> {
        match self.inner.analytic_propagate(&epoch.inner) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to propagate rock: {}", e))),
        }
    }

    fn analytic_at(&self, epoch: PyRef<PyTime>) -> PyResult<PySpaceRock> {
        match self.inner.analytic_at(&epoch.inner) {
            Ok(rock) => Ok(PySpaceRock { inner: rock }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create propagated rock: {}", e))),
        }
    }

    /// Change the reference plane of the SpaceRock.
    ///
    /// # Arguments
    ///
    /// * `reference_plane` - The new reference plane.
    fn change_reference_plane(&mut self, reference_plane: &str) -> PyResult<()> {
        match self.inner.change_reference_plane(reference_plane) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to change reference plane: {}", e))),
        }
    }

    fn to_ssb(&mut self, kernel: &PySpiceKernel) -> PyResult<()> {
        let ker = &kernel.inner;
        match self.inner.to_ssb(&ker) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to change origin to SSB: {}", e))),
        }
    }

    fn to_helio(&mut self, kernel: &PySpiceKernel) -> PyResult<()> {
        let ker = &kernel.inner;
        match self.inner.to_helio(&ker) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to change origin to Heliocenter: {}", e))),
        }
    }

    
    fn __repr__(&self) -> String {
        // build a string representation of the object
        let s = format!("SpaceRock: {}\nposition: {:?}\nvelocity: {:?}\nepoch: {}\nreference_plane: {}\norigin: {}\n", 
            self.inner.name, self.inner.position, self.inner.velocity, self.inner.epoch, self.inner.reference_plane, self.inner.origin);
        s
    }

    #[getter]
    fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    fn reference_plane(&self) -> String {
        self.inner.reference_plane.to_string()
    }

    #[getter]
    fn origin(&self) -> PyOrigin {  
        PyOrigin { inner: self.inner.origin.clone() }
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.to_string()
    }

    #[setter]
    fn set_name(&mut self, name: &str) -> PyResult<()> {
        self.inner.name = name.to_string().into();
        Ok(())
    }

    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        (self.inner.position.x, self.inner.position.y, self.inner.position.z)
    }

    #[setter]
    fn set_position(&mut self, pos: (f64, f64, f64)) -> PyResult<()> {
        self.inner.position = Vector3::new(pos.0, pos.1, pos.2);
        Ok(())
    }

    #[getter]
    fn velocity(&self) -> (f64, f64, f64) {
        (self.inner.velocity.x, self.inner.velocity.y, self.inner.velocity.z)
    }

    #[setter]
    fn set_velocity(&mut self, vel: (f64, f64, f64)) -> PyResult<()> {
        self.inner.velocity = Vector3::new(vel.0, vel.1, vel.2);
        Ok(())
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.position.x
    }

    #[setter]
    fn set_x(&mut self, x: f64) -> PyResult<()> {
        self.inner.position.x = x;
        Ok(())
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.position.y
    }

    #[setter]
    fn set_y(&mut self, y: f64) -> PyResult<()> {
        self.inner.position.y = y;
        Ok(())
    }

    #[getter]
    fn z(&self) -> f64 {
        self.inner.position.z
    }

    #[setter]
    fn set_z(&mut self, z: f64) -> PyResult<()> {
        self.inner.position.z = z;
        Ok(())
    }

    #[getter]
    fn vx(&self) -> f64 {
        self.inner.velocity.x
    }

    #[setter]
    fn set_vx(&mut self, vx: f64) -> PyResult<()> {
        self.inner.velocity.x = vx;
        Ok(())
    }

    #[getter]
    fn vy(&self) -> f64 {
        self.inner.velocity.y
    }

    #[setter]
    fn set_vy(&mut self, vy: f64) -> PyResult<()> {
        self.inner.velocity.y = vy;
        Ok(())
    }

    #[getter]
    fn vz(&self) -> f64 {
        self.inner.velocity.z
    }

    #[setter]
    fn set_vz(&mut self, vz: f64) -> PyResult<()> {
        self.inner.velocity.z = vz;
        Ok(())
    }

    fn set_mass(&mut self, mass: f64) -> PyResult<()> {
        self.inner.set_mass(mass);
        Ok(())
    }

    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass()
    }

    fn set_absolute_magnitude(&mut self, absolute_magnitude: f64) -> PyResult<()> {
        self.inner.set_absolute_magnitude(absolute_magnitude);
        Ok(())
    }

    #[getter]
    fn absolute_magnitude(&self) -> f64 {
        self.inner.absolute_magnitude()
    }

    fn set_gslope(&mut self, gslope: f64) -> PyResult<()> {
        self.inner.set_gslope(gslope);
        Ok(())
    }

    #[getter]
    fn gslope(&self) -> f64 {
        self.inner.gslope()
    }

    #[getter]
    fn evec(&self) -> (f64, f64, f64) {
        let e = self.inner.evec();
        (e.x, e.y, e.z)
    }

    #[getter]
    fn mu(&self) -> f64 {
        self.inner.origin.mu()
    }

    pub fn a(&self) -> f64 {
        self.inner.a()
    }

    pub fn q(&self) -> f64 {
        self.inner.q()
    }

    pub fn p(&self) -> f64 {
        self.inner.p()
    }

    pub fn e(&self) -> f64 {
        self.inner.e()
    }

    pub fn inc(&self) -> f64 {
        self.inner.inc()
    }

    pub fn arg(&self) -> f64 {
        self.inner.arg()
    }

    pub fn node(&self) -> f64 {
        self.inner.node()
    }

    pub fn true_anomaly(&self) -> f64 {
        self.inner.true_anomaly()
    }

    pub fn mean_anomaly(&self) -> f64 {
        self.inner.mean_anomaly()
    }

    pub fn conic_anomaly(&self) -> f64 {
        self.inner.conic_anomaly()
    }



}

