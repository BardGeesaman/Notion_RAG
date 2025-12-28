"""Bayesian dose-response fitting utilities (PyMC).

Implements a 4-parameter logistic (4PL) model:
  y = bottom + (top - bottom) / (1 + (x/ec50)^(-hill))
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

# Import PriorConfig from analysis layer (no circular import)
from amprenta_rag.analysis.models import PriorConfig


def _require_pymc():
    try:
        import pymc as pm  # type: ignore
        import arviz as az  # type: ignore
        import numpy as np  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "PyMC/ArviZ is required for Bayesian dose-response. "
            "Install with: pymc>=5.0 arviz>=0.16 pytensor>=2.14"
        ) from e
    return pm, az, np


def fit_bayesian_dose_response(
    concentrations: List[float],
    responses: List[float],
    prior_ec50: Optional[float] = None,
    likelihood: str = "normal",  # or "student_t"
    prior_config: Optional["PriorConfig"] = None,
) -> Dict[str, Any]:
    """
    Fit a 4-parameter logistic dose-response model with PyMC.

    Returns:
        {
          "ec50_mean": float,
          "ec50_ci": (low, high),
          "hill_slope": float,
          "trace": arviz.InferenceData,
          "posterior_predictive": arviz.InferenceData,
        }
    """

    pm, az, np = _require_pymc()

    if len(concentrations) != len(responses):
        raise ValueError("concentrations and responses must have the same length")
    if len(concentrations) < 6:
        raise ValueError("Need at least 6 points for Bayesian dose-response fit")

    x = np.asarray(concentrations, dtype=float)
    y = np.asarray(responses, dtype=float)

    if np.any(x <= 0):
        raise ValueError("All concentrations must be > 0 for 4PL model")

    y_min = float(np.nanmin(y))
    y_max = float(np.nanmax(y))
    y_range = max(1e-8, y_max - y_min)

    # Prior EC50 default: median concentration.
    if prior_ec50 is None:
        prior_ec50 = float(np.median(x))

    like = (likelihood or "normal").strip().lower()
    if like not in {"normal", "student_t"}:
        raise ValueError("likelihood must be 'normal' or 'student_t'")

    with pm.Model():
        # Configure priors based on prior_config or use defaults
        if prior_config:
            # Bottom prior
            bottom_mean = prior_config.bottom_prior_mean if prior_config.bottom_prior_mean is not None else y_min
            bottom = pm.Normal("bottom", mu=bottom_mean, sigma=2.0 * y_range)
            
            # Top prior  
            top_mean = prior_config.top_prior_mean if prior_config.top_prior_mean is not None else y_max
            top = pm.Normal("top", mu=top_mean, sigma=2.0 * y_range)
            
            # EC50 prior (log-scale)
            ec50_mean = prior_config.ec50_prior_mean if prior_config.ec50_prior_mean is not None else np.log(max(1e-12, prior_ec50))
            ec50 = pm.LogNormal("ec50", mu=ec50_mean, sigma=prior_config.ec50_prior_sd)
            
            # Hill slope prior
            hill_slope = pm.Normal("hill_slope", mu=prior_config.hill_prior_mean, sigma=prior_config.hill_prior_sd)
        else:
            # Default priors (original behavior)
            bottom = pm.Normal("bottom", mu=y_min, sigma=2.0 * y_range)
            top = pm.Normal("top", mu=y_max, sigma=2.0 * y_range)
            ec50 = pm.LogNormal("ec50", mu=np.log(max(1e-12, prior_ec50)), sigma=1.0)
            hill_slope = pm.Normal("hill_slope", mu=1.0, sigma=2.0)

        mu = bottom + (top - bottom) / (1.0 + (x / ec50) ** (-hill_slope))

        if like == "normal":
            sigma = pm.HalfNormal("sigma", sigma=y_range)
            pm.Normal("y", mu=mu, sigma=sigma, observed=y)
        else:
            sigma = pm.HalfNormal("sigma", sigma=y_range)
            nu = pm.Exponential("nu", lam=1 / 10)
            pm.StudentT("y", nu=nu, mu=mu, sigma=sigma, observed=y)

        trace = pm.sample(
            draws=300,
            tune=300,
            chains=2,
            target_accept=0.9,
            progressbar=False,
            random_seed=42,
        )

        ppc = pm.sample_posterior_predictive(trace, var_names=["y"], progressbar=False)

    # Summaries
    ec50_samples = trace.posterior["ec50"].values.reshape(-1)
    hill_samples = trace.posterior["hill_slope"].values.reshape(-1)

    ec50_mean = float(np.mean(ec50_samples))
    ec50_ci = (float(np.quantile(ec50_samples, 0.025)), float(np.quantile(ec50_samples, 0.975)))
    hill_mean = float(np.mean(hill_samples))

    return {
        "ec50_mean": ec50_mean,
        "ec50_ci": ec50_ci,
        "hill_slope": hill_mean,
        "trace": trace,
        "posterior_predictive": ppc,
    }


__all__ = ["fit_bayesian_dose_response"]


