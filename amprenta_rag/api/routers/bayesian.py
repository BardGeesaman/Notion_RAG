"""Bayesian optimization API endpoints."""
from __future__ import annotations

from typing import Dict, List, Optional

from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective

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


__all__ = ["router"]

