use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::assist::forces::{
    assist_default_forces, EarthHarmonics, Force, GrEih, GrPotential, GrSimple, NewtonianGravity, NonGravitational, SolarJ2,
};
use spacerocks::assist::EphemerisConstants;

use crate::py_spice::spicekernel::PySpiceKernel;

fn constants(kernel: Option<PyRef<PySpiceKernel>>) -> EphemerisConstants {
    match kernel {
        Some(k) => EphemerisConstants::from_kernel(&k.inner),
        None => EphemerisConstants::DE440,
    }
}

/// A force model for `SpiceSimulation` (the models of ASSIST). Constants (J2, radii, speed of
/// light) are read from `kernel`'s planetary ephemeris when given, DE440 values otherwise.
#[pyclass]
#[pyo3(name = "Force")]
pub struct PyForce {
    pub inner: Box<dyn Force + Send + Sync>,
    pub name: String,
}

impl PyForce {
    fn new(inner: Box<dyn Force + Send + Sync>, name: &str) -> Self {
        PyForce { inner, name: name.to_string() }
    }
}

#[pymethods]
impl PyForce {
    /// Newtonian point-mass gravity of the simulation's SPICE perturbers.
    #[classmethod]
    pub fn newtonian_gravity(_cls: Py<PyType>) -> Self {
        PyForce::new(Box::new(NewtonianGravity), "newtonian_gravity")
    }

    /// Earth's J2, J3, J4 (pole fixed at the J2000 pole, as in ASSIST and Horizons).
    #[classmethod]
    #[pyo3(signature = (kernel = None))]
    pub fn earth_harmonics(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>) -> Self {
        PyForce::new(Box::new(EarthHarmonics::new(&constants(kernel))), "earth_harmonics")
    }

    /// The Sun's J2.
    #[classmethod]
    #[pyo3(signature = (kernel = None))]
    pub fn solar_j2(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>) -> Self {
        PyForce::new(Box::new(SolarJ2::new(&constants(kernel))), "solar_j2")
    }

    /// Einstein-Infeld-Hoffmann relativity, with the first `sources` of Sun, Mercury, Venus,
    /// Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto as sources (ASSIST default: 1).
    #[classmethod]
    #[pyo3(signature = (kernel = None, sources = 1))]
    pub fn gr_eih(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>, sources: usize) -> Self {
        PyForce::new(Box::new(GrEih::new(&constants(kernel)).with_sources(sources)), "gr_eih")
    }

    /// The Sun's one-body post-Newtonian correction.
    #[classmethod]
    #[pyo3(signature = (kernel = None))]
    pub fn gr_simple(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>) -> Self {
        PyForce::new(Box::new(GrSimple::new(&constants(kernel))), "gr_simple")
    }

    /// A velocity-independent approximation to the Sun's relativistic correction.
    #[classmethod]
    #[pyo3(signature = (kernel = None))]
    pub fn gr_potential(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>) -> Self {
        PyForce::new(Box::new(GrPotential::new(&constants(kernel))), "gr_potential")
    }

    /// Marsden-model non-gravitational accelerations A_k g(r), with
    /// g(r) = alpha (r/r0)^-m (1 + (r/r0)^n)^-k. The defaults are JPL's asteroid (1/r^2) law;
    /// each particle's A1, A2, A3 come from `SpaceRock.set_nongrav`.
    #[classmethod]
    #[pyo3(signature = (alpha = 1.0, r0 = 1.0, m = 2.0, n = 5.093, k = 0.0))]
    pub fn nongravitational(_cls: Py<PyType>, alpha: f64, r0: f64, m: f64, n: f64, k: f64) -> Self {
        PyForce::new(Box::new(NonGravitational { alpha, r0, m, n, k }), "nongravitational")
    }

    /// Non-gravitational accelerations with the standard comet (water-ice) g(r).
    #[classmethod]
    pub fn nongravitational_comet(_cls: Py<PyType>) -> Self {
        PyForce::new(Box::new(NonGravitational::comet()), "nongravitational_comet")
    }

    /// ASSIST's default force model (the one `SpiceSimulation.horizons` uses), as a list.
    #[classmethod]
    #[pyo3(signature = (kernel = None))]
    pub fn assist_defaults(_cls: Py<PyType>, kernel: Option<PyRef<PySpiceKernel>>) -> Vec<PyForce> {
        let names = ["nongravitational", "earth_harmonics", "solar_j2", "gr_eih", "newtonian_gravity"];
        assist_default_forces(&constants(kernel))
            .into_iter()
            .zip(names)
            .map(|(f, n)| PyForce::new(f, n))
            .collect()
    }

    fn __repr__(&self) -> String {
        format!("Force.{}", self.name)
    }
}
