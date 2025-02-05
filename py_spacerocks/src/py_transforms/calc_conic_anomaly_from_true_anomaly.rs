use pyo3::prelude::*;

use spacerocks::transforms::calc_conic_anomaly_from_true_anomaly;

#[pyfunction]
#[pyo3(name = "calc_conic_anomaly_from_true_anomaly")]
/// Calculate the conic anomaly from the true anomaly.
///
/// # Arguments
/// * `e` Eccentricity of the conic section.
/// * `true_anomaly` True anomaly in radians.
///
/// # Returns
/// * `conic_anomaly` Conic anomaly in radians.
pub fn calc_conic_anomaly_from_true_anomaly_py(e: f64, true_anomaly: f64) -> PyResult<f64> {
    match calc_conic_anomaly_from_true_anomaly(e, true_anomaly) {
        Ok(result) => Ok(result),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e))),
    }
}