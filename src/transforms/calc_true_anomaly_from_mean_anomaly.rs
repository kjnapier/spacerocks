use crate::transforms::{calc_conic_anomaly_from_mean_anomaly, calc_true_anomaly_from_conic_anomaly};
use crate::errors::OrbitError;

/// Calculate the true anomaly from the mean anomaly
pub fn calc_true_anomaly_from_mean_anomaly(e: f64, mean_anomaly: f64) -> Result<f64, OrbitError> {
    let conic_anomaly = calc_conic_anomaly_from_mean_anomaly(e, mean_anomaly)?;
    let true_anomaly = calc_true_anomaly_from_conic_anomaly(e, conic_anomaly)?; 
    Ok(true_anomaly)
}