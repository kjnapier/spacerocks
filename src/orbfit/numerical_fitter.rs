pub struct NumericalOrbitFitter {
    pub observations: Vec<Observation>,
    pub simulation: Simulation,
}

impl OrbitFitter for NumericalOrbitFitter {
    type Jacobian = Vec<Vec<f64>>;

    fn residuals(&self, params: &OrbitParameters) -> Vec<f64> {
        self.observations.iter().map(|obs| {
            let predicted = self.simulation.simulate(params, obs.time);
            // For illustration, use the difference in the x-component.
            predicted.0 - obs.position.0
        }).collect()
    }

    fn jacobian(&self, params: &OrbitParameters) -> Self::Jacobian {
        let base_residuals = self.residuals(params);
        let n_params = 6; // Assuming six orbit parameters.
        let n_res = base_residuals.len();
        let mut jacobian = vec![vec![0.0; n_params]; n_res];

        for i in 0..n_params {
            let mut params_perturbed = params.clone();
            perturb_parameter(&mut params_perturbed, i, self.epsilon);
            let perturbed_residuals = self.residuals(&params_perturbed);
            for j in 0..n_res {
                jacobian[j][i] = (perturbed_residuals[j] - base_residuals[j]) / self.epsilon;
            }
        }

        jacobian
    }
}
