
1. Have a list of Detections
2. Posit a guess for the orbit
    - Gauss' Method, for example
3. Choose a model for the orbit
    - Keplerian, N-body, etc.
4. Use the model to calculate residuals
5. Minimize the residuals
    - Levenberg-Marquardt, for example


THE DATA AND ALGORITHMS NEED TO REMAIN SEPARATE!!! 


Model is a trait that can be implemented for any orbit model. It should have the following methods:
    - residuals
    - jacobian
    - residuals_and_jacobian
    - chi_squared
    
    The Model needs to be instantiated with either nothing (for a Keplerian model) or with a Simulation object (for an N-body model).
    - KeplerianModel
    - NBodyModel(Simulation)

    The behavior of each model type is defines inside of the implementation of the Model trait.

Then there will be a Minimizer trait that can be implemented for any minimization algorithm. 

    Different minimizers will need to be instantiates with different parameters. For example, the Levenberg-Marquardt minimizer will need to be instantiated with a damping factor, while the Nelder-Mead minimizer will not.

    minimizer = LevenbergMarquardt(tol=1e-6, lambda=1e-3, ...)

    It should just take a model and a list of detections and return the best-fit model.
    minimizer.fit(model, detections, initial_guess)

    It should also implement a `step` method that can be used to take a single step in the minimization process.
    minimizer.step(model, detections)

    The behavior of each minimization algorithm is defined inside of the implementation of the Minimizer trait.

