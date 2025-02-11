/// Trait defining the interface for orbital models that can be fitted to data
pub trait Model {
    type Jacobian;

    /// Compute residuals between model predictions and observations
    ///
    /// # Arguments
    /// * `params` - Current orbital parameters
    ///
    /// # Returns
    /// Vector of residuals between model and observations
    fn residuals(&self, params: &OrbitParameters) -> Vec<f64>;

    /// Compute the Jacobian matrix for the given orbit parameters.
    ///
    /// # Arguments
    /// * `params` - Current orbital parameters
    ///
    /// # Returns
    /// Jacobian matrix of partial derivatives
    fn jacobian(&self, params: &OrbitParameters) -> Self::Jacobian;
}