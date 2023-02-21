from typing import List, Callable, Optional
#from detection import Detection

import numpy as np

def detprob_logit(m, params: List[float]):
    m50, eta0, sigma = params
    logit = (eta0 / 2) * (1 - np.tanh((m - m50) / sigma))
    return logit

def detprob_double(m, params: List[float]):
    m25, c, k1, k2 = params
    logit = c / ((1 + np.exp((m - m25)/k1)) * (1 + np.exp((m - m25)/k2)))
    return logit

def logit_prior(theta: List[float]):
    m50, eta0, sigma = theta
    if not (0 <= eta0 <= 1):
        return np.inf
    return 0

def double_prior(theta: List[float]):
    m50, eta0, k1, k2 = theta
    if not (0 <= eta0 <= 1):
        return np.inf
    return 0

def ossos_prior(theta: List[float]):
    m0, eta0, c, sigma = params = theta
    if not (0 <= eta0 <= 1):
        return np.inf
    return 0

def logL(theta: List[float], magnitudes: List[float], weights: List[float], form: Callable):
    p = form(magnitudes, theta)
    return np.sum(weights * np.log(p)) + np.sum((1 - weights) * np.log(1 - p))
  
def minuslogL(theta: List[float], magnitudes: List[float], weights: List[float], form: Callable):
    return -1 * logL(theta, magnitudes, weights, form)

def log_probability(theta: List[float], magnitudes: List[float], weights: List[float], form: Callable, log_prior: Callable):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    val = lp + logL(theta, magnitudes, weights, form)
    if np.isnan(val):
        return -np.inf
    return val

def truncated_efficiency(m, form, cut):
    m = np.atleast_1d(m)
    etas = np.zeros_like(m)
    keep = m < cut
    etas[keep] = form(m[keep])
    etas[~keep] = 0
    return etas

efficiency_forms = {'logit': [detprob_logit, 
                              logit_prior, 
                              ((25, 27), (0, 1), (0, 100)), 
                              (26, 1, 1), 
                              ['m_{50}', '\eta_0', '\sigma']],
                    'double': [detprob_double, 
                               double_prior, 
                               ((25, 27), (0, 1), (0, 100), (0, 100)), 
                               (26, 1, 1, 1), 
                               ['m_{25}', '\eta_0', 'k_1', 'k_2']]}

import functools
from scipy.optimize import minimize
import pandas as pd

class EfficiencyFunction:

    def __init__(self, form: Callable, labels: List[str], prior: Optional[Callable] = None):
        self.form = form
        self.labels = labels
        self.prior = prior

    def fit(self, magnitude: List[float], weight: List[float], bounds: List[tuple], theta0: List[float]):
        res = minimize(minuslogL, 
                       x0=theta0, 
                       method='Powell', 
                       args=(magnitude, weight, self.form), 
                       tol=1e-8,
                       bounds=bounds)
        return functools.partial(self.form, params=res['x']), res, {l: x for l, x in zip(logit.labels, res['x'])}

logit = EfficiencyFunction(detprob_logit, ['m_{50}', '\eta_0', '\sigma'])
double = EfficiencyFunction(detprob_double, ['m_{25}', '\eta_0', 'k_1', 'k_2'])