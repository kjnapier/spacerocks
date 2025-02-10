pub trait Model {
    type Jacobian;

    /// Compute the residual vector for the given orbit parameters.
    fn residuals(&self, params: &OrbitParameters) -> Vec<f64>;

    /// Compute the Jacobian matrix for the given orbit parameters.
    fn jacobian(&self, params: &OrbitParameters) -> Self::Jacobian;
}