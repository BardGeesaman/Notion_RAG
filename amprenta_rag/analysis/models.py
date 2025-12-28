"""Analysis layer models and data structures."""
from __future__ import annotations

from typing import Optional
from pydantic import BaseModel


class PriorConfig(BaseModel):
    """Configuration for Bayesian model priors."""
    
    ec50_prior_mean: Optional[float] = None  # Log-scale, None = auto from data
    ec50_prior_sd: float = 1.0
    hill_prior_mean: float = 1.0
    hill_prior_sd: float = 2.0
    bottom_prior_mean: Optional[float] = None  # None = auto from data
    top_prior_mean: Optional[float] = None
