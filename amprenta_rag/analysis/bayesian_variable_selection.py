"""Bayesian variable selection with Horseshoe prior for biomarker discovery."""
from __future__ import annotations

from typing import Any, Dict, List
import numpy as np


def _require_pymc():
    try:
        import pymc as pm
        import arviz as az
    except ImportError as e:
        raise ImportError("PyMC/ArviZ required. Install: pymc>=5.0 arviz>=0.16") from e
    return pm, az


def compute_posterior_inclusion_probabilities(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: List[str],
    pip_threshold: float = 0.5,
    n_samples: int = 1000,
) -> Dict[str, Any]:
    """
    Compute Posterior Inclusion Probabilities using Horseshoe prior.
    
    Args:
        X: Feature matrix (n_samples, n_features)
        y: Response vector (n_samples,)
        feature_names: Names for each feature column
        pip_threshold: Threshold for "included" classification (default 0.5)
        n_samples: MCMC samples (default 1000)
    
    Returns:
        {
            "pips": {feature_name: probability},
            "included_features": [names with PIP > threshold],
            "coefficients": {feature_name: posterior_mean},
            "summary": arviz summary DataFrame
        }
    """
    pm, az = _require_pymc()
    
    n_obs, n_feat = X.shape
    if len(feature_names) != n_feat:
        raise ValueError(f"feature_names length ({len(feature_names)}) != X columns ({n_feat})")
    
    # Standardize
    X_std = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)
    y_std = (y - y.mean()) / (y.std() + 1e-8)
    
    with pm.Model():
        # Horseshoe prior (manual implementation for PyMC 5.0+ compatibility)
        tau = pm.HalfCauchy("tau", beta=1.0)  # Global shrinkage
        lambdas = pm.HalfCauchy("lambdas", beta=1.0, shape=n_feat)  # Local shrinkage
        
        beta = pm.Normal("beta", mu=0, sigma=tau * lambdas, shape=n_feat)
        
        mu = pm.math.dot(X_std, beta)
        sigma = pm.HalfNormal("sigma", sigma=1.0)
        pm.Normal("y", mu=mu, sigma=sigma, observed=y_std)
        
        trace = pm.sample(draws=n_samples, tune=500, chains=2, 
                         progressbar=False, random_seed=42)
    
    # Compute PIPs: P(|beta| > threshold)
    beta_samples = trace.posterior["beta"].values.reshape(-1, n_feat)
    pips = {}
    coeffs = {}
    
    for i, name in enumerate(feature_names):
        samples = beta_samples[:, i]
        # PIP = proportion of samples where |beta| > small threshold
        pips[name] = float(np.mean(np.abs(samples) > 0.01))
        coeffs[name] = float(np.mean(samples))
    
    # Sort by PIP
    pips = dict(sorted(pips.items(), key=lambda x: x[1], reverse=True))
    included = [k for k, v in pips.items() if v > pip_threshold]
    
    return {
        "pips": pips,
        "included_features": included,
        "coefficients": coeffs,
        "n_included": len(included),
        "pip_threshold": pip_threshold,
    }


__all__ = ["compute_posterior_inclusion_probabilities"]
