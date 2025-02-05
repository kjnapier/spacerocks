use pyo3::prelude::*;

use spacerocks::transforms::calc_true_anomaly_from_conic_anomaly;

#[pyfunction]
#[pyo3(name = "calc_true_anomaly_from_conic_anomaly")]
/// Calculate the true anomaly from the conic anomaly.
///
/// # Arguments
/// * `e` Eccentricity of the conic section.
/// * `conic_anomaly` Conic anomaly in radians.
///
/// # Returns
/// * `true_anomaly` True anomaly in radians.
pub fn calc_true_anomaly_from_conic_anomaly_py(e: f64, conic_anomaly: f64) -> PyResult<f64> {
    match calc_true_anomaly_from_conic_anomaly(e, conic_anomaly) {
        Ok(result) => Ok(result),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e))),
    }
}