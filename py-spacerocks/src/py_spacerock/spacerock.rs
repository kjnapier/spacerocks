use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::spacerock::SpaceRock;
use spacerocks::spacerock::CoordinateFrame;

use spacerocks::KeplerOrbit;
use spacerocks::Properties;
use spacerocks::constants::MU_BARY;
use spacerocks::transforms::calc_xyz_from_kepM;

use nalgebra::Vector3;

use crate::py_time::time::PyTime;
use crate::py_observing::detection::PyDetection;
use crate::py_observing::observer::PyObserver;
use crate::py_spacerock::origin::PyOrigin;


#[pyclass]
#[pyo3(name = "SpaceRock")]
pub struct PySpaceRock {
    pub inner: SpaceRock,
}

#[pymethods]
impl PySpaceRock {

    #[classmethod]
    // #[pyo3(signature = (name, epoch, frame="J2000", origin=PyOrigin::ssb()))]
    fn from_horizons(_cls: &PyType, name: &str, epoch: PyRef<PyTime>, frame: &str, origin: PyRef<PyOrigin>) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let rock = SpaceRock::from_horizons(name, ep, &frame, or);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Horizons for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    // #[pyo3(signature = (name, epoch, frame="J2000", origin=PyOrigin::ssb()))]
    fn from_spice(_cls: &PyType, name: &str, epoch: &PyTime, frame: &str, origin: &PyOrigin) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let rock = SpaceRock::from_spice(name, &ep, &frame, &or);
        Ok(PySpaceRock { inner: rock })
    }

    #[classmethod]
    fn random(_cls: &PyType) -> Self {
        PySpaceRock { inner: SpaceRock::random() }
    }

    #[classmethod]
    fn from_xyz(_cls: &PyType, name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: PyRef<PyTime>, frame: &str, origin: PyRef<PyOrigin>) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, ep.clone(), &frame, or);
        Ok(PySpaceRock { inner: rock })
    }

    #[classmethod]
    fn from_spherical(_cls: &PyType, name: &str, phi: f64, theta: f64, r: f64, vr: f64, vo: f64, psi: f64, epoch: PyRef<PyTime>, frame: &str, origin: PyRef<PyOrigin>) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let rock = SpaceRock::from_spherical(name, phi, theta, r, vr, vo, psi, ep.clone(), &frame, or);
        Ok(PySpaceRock { inner: rock })
    }

    #[classmethod]
    // #[pyo3(signature = (name, a, e, inc, arg, node, M, epoch, frame="ECLIPJ2000", origin=PyOrigin::ssb()))]
    fn from_kepler(_cls: &PyType, name: &str, a: f64, e: f64, inc: f64, arg: f64, node: f64, M: f64, epoch: PyRef<PyTime>, frame: &str, origin: PyRef<PyOrigin>) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;
        let state = calc_xyz_from_kepM(a, e, inc, arg, node, M);
        let acceleration = Vector3::new(0.0, 0.0, 0.0);

        // let frame = CoordinateFrame::from_str(frame).unwrap();
        let frame = CoordinateFrame::from_str(frame).unwrap();

        let rock = SpaceRock {
            name: name.to_string().into(),
            position: state.position,
            velocity: state.velocity,
            acceleration: acceleration,
            epoch: ep.clone(),
            frame: frame,
            origin: or.clone(),
            orbit: None,
            mass: 0.0,
            properties: None,
        };
        Ok(PySpaceRock { inner: rock })
    }

    fn set_absolute_magnitude(&mut self, H: f64) {
        let properties = Properties {
            H: Some(H),
            Gslope: Some(0.15),
            radius: None,
            albedo: None,
        };
        self.inner.properties = Some(properties);
    }
    
    fn __repr__(&self) -> String {
        format!("SpaceRock: {}", self.inner.name)
    }

    fn analytic_propagate(&mut self, t: &PyTime) {
        self.inner.analytic_propagate(&t.inner);
    }

    fn at(&mut self, t: &PyTime) {
        self.inner.at(&t.inner);
    }

    fn observe(&mut self, observer: &PyObserver) -> PyResult<PyDetection> {
        if observer.inner.frame != CoordinateFrame::J2000 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks.")));
        }

        match self.inner.observe(&observer.inner) {
            Ok(obs) => Ok(PyDetection { inner: obs }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks. Error: {}", e))),
        }
    }

    fn change_frame(&mut self, frame: &str) -> PyResult<()> {
        self.inner.change_frame(frame);
        Ok(())
    }

    fn calculate_orbit(&mut self) {
        self.inner.calculate_orbit();
    }

    #[getter]
    fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    fn frame(&self) -> String {
        self.inner.frame.to_string()
    }

    #[getter]
    fn origin(&self) -> String {
        self.inner.origin.to_string()   
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    #[getter]
    fn e(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.e),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn a(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.a),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn inc(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.inc),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn node(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.node),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn arg(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.arg),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn varpi(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.varpi()),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn name(&self) -> String {
        (*self.inner.name).clone()
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

    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass
    }

    #[setter]
    fn set_mass(&mut self, mass: f64) -> PyResult<()> {
        self.inner.mass = mass;
        Ok(())
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

    #[getter]
    fn mean_anomaly(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.M()),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn mean_longitude(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok((orbit.M() + orbit.varpi()) % (2.0 * std::f64::consts::PI)),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

    #[getter]
    fn t_peri(&self) -> PyResult<PyTime> {
        let t0 = &self.inner.epoch;
        let mean_motion = (self.a()?.powi(3) / MU_BARY).sqrt();
        let dt = self.mean_anomaly()? * mean_motion;
        Ok(PyTime { inner: t0 + dt })
    }

    #[getter]
    fn q(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.q()),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated. Use the calculate_orbit() method to calculate the orbit.")),
        }
    }

}

