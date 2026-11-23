"""Small least-squares helpers.

The reports quote fit uncertainties as standard errors of the parameters.
Two things are done more carefully here:
  * derived quantities (e.g. a threshold -b/s) are either fitted directly or
    propagated with the full covariance matrix, since slope and intercept of a
    line far from x = 0 are strongly anti-correlated;
  * for unweighted fits the parameter covariance is scaled by the residual
    variance, which is what "standard error of the fit" means.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit


@dataclass
class FitResult:
    params: np.ndarray
    cov: np.ndarray
    residuals: np.ndarray
    chi2: float
    dof: int

    @property
    def errors(self):
        return np.sqrt(np.diag(self.cov))

    @property
    def reduced_chi2(self):
        return self.chi2 / self.dof if self.dof > 0 else np.nan


def fit(model, x, y, p0, sigma=None):
    """Wrapper around curve_fit that also returns residuals and chi^2.

    With sigma given, the uncertainties are taken as absolute (chi^2 is then
    a real goodness-of-fit measure). Without sigma, the covariance is scaled
    by the residual variance.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    absolute = sigma is not None
    params, cov = curve_fit(model, x, y, p0=p0, sigma=sigma, absolute_sigma=absolute)
    res = y - model(x, *params)
    weights = 1.0 if sigma is None else 1.0 / np.asarray(sigma, float) ** 2
    chi2 = float(np.sum(weights * res**2))
    return FitResult(params, cov, res, chi2, len(x) - len(params))


def line(x, m, c):
    return m * x + c


def hinge(x, slope, x_th):
    """Straight line written through its zero crossing: slope * (x - x_th)."""
    return slope * (x - x_th)


def propagate_ratio(m, c, cov):
    """Zero crossing x0 = -c/m of a line and its error, with and without
    the covariance term. Returns (x0, sigma_full, sigma_no_cov)."""
    x0 = -c / m
    grad = np.array([c / m**2, -1.0 / m])  # d x0 / d(m, c)
    var_full = grad @ cov @ grad
    var_diag = grad**2 @ np.diag(cov)
    return x0, np.sqrt(var_full), np.sqrt(var_diag)


def residual_bootstrap(model, x, y, p0, n_boot=4000, rng=None):
    """Refit after resampling residuals; returns the array of parameter sets."""
    rng = np.random.default_rng(rng)
    base = fit(model, x, y, p0)
    y_hat = model(np.asarray(x, float), *base.params)
    res = base.residuals - base.residuals.mean()
    out = np.empty((n_boot, len(p0)))
    for i in range(n_boot):
        y_star = y_hat + rng.choice(res, size=len(res), replace=True)
        out[i] = curve_fit(model, x, y_star, p0=base.params)[0]
    return out
