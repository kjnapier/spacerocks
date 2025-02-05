use pyo3::prelude::*;

use spacerocks::transforms::calc_true_anomaly_from_mean_anomaly;

#[pyfunction]
#[pyo3(name = "calc_true_anomaly_from_mean_anomaly")]
/// Calculate the true anomaly from the mean anomaly.
///
/// # Arguments
/// * `e` Eccentricity of the conic section.
/// * `mean_anomaly` Mean anomaly in radians.
///
/// # Returns
/// * `true_anomaly` True anomaly in radians.
pub fn calc_true_anomaly_from_mean_anomaly_py(e: f64, mean_anomaly: f64) -> PyResult<f64> {
    match calc_true_anomaly_from_mean_anomaly(e, mean_anomaly) {
        Ok(result) => Ok(result),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e))),
    }
}