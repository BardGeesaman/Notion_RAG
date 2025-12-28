"""Bayesian optimization API endpoints."""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
from amprenta_rag.analysis.bayesian_dose_response import fit_bayesian_dose_response
from amprenta_rag.api.schemas import PriorConfig

router = APIRouter(prefix="/bayesian-optimization", tags=["bayesian-optimization"])


@router.post("/recommend-multi-objective")
def recommend_multi_objective_endpoint(
    tested_compounds: List[Dict],
    candidate_pool: List[Dict],
    objectives: List[str],
    objective_directions: List[str],
    batch_size: int = 10,
    ref_point: Optional[List[float]] = None,
) -> Dict:
    """Multi-objective Bayesian optimization recommendations using qNEHVI.
    
    Args:
        tested_compounds: List of tested compounds with features and objective values
        candidate_pool: List of candidate compounds with id and features
        objectives: List of objective names to optimize
        objective_directions: List of "maximize" or "minimize" for each objective
        batch_size: Number of recommendations to return
        ref_point: Optional reference point for hypervolume calculation
        
    Returns:
        Dict with recommendations and Pareto front information
        
    Raises:
        HTTPException: 400 for validation errors, 500 for computation errors
    """
    try:
        return recommend_multi_objective(
            tested_compounds=tested_compounds,
            candidate_pool=candidate_pool,
            objectives=objectives,
            objective_directions=objective_directions,
            batch_size=batch_size,
            ref_point=ref_point,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Optimization failed: {e}")


@router.post("/dose-response")
def fit_dose_response_endpoint(
    concentrations: List[float],
    responses: List[float],
    prior_config: Optional[PriorConfig] = None,
    likelihood: str = "normal",
) -> Dict:
    """Fit Bayesian dose-response model with configurable priors.
    
    Args:
        concentrations: List of concentration values (must be > 0)
        responses: List of response values
        prior_config: Optional prior configuration for Bayesian fitting
        likelihood: "normal" or "student_t" likelihood function
        
    Returns:
        Dictionary with fitting results including EC50 estimates and credible intervals
    """
    try:
        if len(concentrations) != len(responses):
            raise HTTPException(
                status_code=400, 
                detail="concentrations and responses must have the same length"
            )
        
        if len(concentrations) < 6:
            raise HTTPException(
                status_code=400,
                detail="Need at least 6 points for Bayesian dose-response fit"
            )
        
        if any(c <= 0 for c in concentrations):
            raise HTTPException(
                status_code=400,
                detail="All concentrations must be > 0 for 4PL model"
            )
        
        if likelihood not in {"normal", "student_t"}:
            raise HTTPException(
                status_code=400,
                detail="likelihood must be 'normal' or 'student_t'"
            )
        
        result = fit_bayesian_dose_response(
            concentrations=concentrations,
            responses=responses,
            likelihood=likelihood,
            prior_config=prior_config,
        )
        
        return result
        
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Dose-response fitting failed: {e}")


@router.post("/variable-selection")
def variable_selection_endpoint(
    features: List[Dict[str, Any]],  # [{"name": str, "values": List[float]}]
    response: List[float],
    pip_threshold: float = 0.5,
    n_samples: int = 1000,
) -> Dict:
    """Bayesian variable selection using Horseshoe prior."""
    try:
        # Convert to numpy arrays
        feature_names = [f["name"] for f in features]
        X = np.column_stack([f["values"] for f in features])
        y = np.array(response)
        
        from amprenta_rag.analysis.bayesian_variable_selection import (
            compute_posterior_inclusion_probabilities
        )
        
        return compute_posterior_inclusion_probabilities(
            X=X, y=y, feature_names=feature_names,
            pip_threshold=pip_threshold, n_samples=n_samples
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Variable selection failed: {e}")


__all__ = ["router"]

