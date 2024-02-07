use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Properties {
    pub H: Option<f64>,
    pub Gslope: Option<f64>,
    pub radius: Option<f64>,
    pub albedo: Option<f64>,
}