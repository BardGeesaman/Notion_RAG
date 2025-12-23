"""Analysis API endpoints (Bayesian inference, optimization, etc.)."""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from amprenta_rag.analysis.bayesian_dose_response import fit_bayesian_dose_response
from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds
from amprenta_rag.api.schemas import BayesianDoseResponseRequest, BayesianDoseResponseResponse


router = APIRouter()

@router.post("/analysis/dose-response/bayesian", response_model=BayesianDoseResponseResponse)
def bayesian_dose_response(payload: BayesianDoseResponseRequest) -> BayesianDoseResponseResponse:
    try:
        out: Dict[str, Any] = fit_bayesian_dose_response(
            concentrations=payload.concentrations,
            responses=payload.responses,
            prior_ec50=payload.prior_ec50,
            likelihood=payload.likelihood,
        )

        diagnostics = None
        if payload.include_diagnostics:
            import base64
            import io

            import arviz as az  # type: ignore
            import matplotlib

            matplotlib.use("Agg")  # Non-interactive backend
            import matplotlib.pyplot as plt

            axes = az.plot_trace(out["trace"], var_names=["ec50", "hill_slope"])
            fig = axes.ravel()[0].figure
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
            plt.close(fig)
            trace_b64 = base64.b64encode(buf.getvalue()).decode("utf-8")

            summary = az.summary(out["trace"], var_names=["ec50", "hill_slope"])
            diagnostics = {
                "trace_plot": f"data:image/png;base64,{trace_b64}",
                "rhat": summary["r_hat"].to_dict(),
                "ess_bulk": summary["ess_bulk"].to_dict(),
            }

        return BayesianDoseResponseResponse(
            ec50_mean=float(out["ec50_mean"]),
            ec50_ci=tuple(out["ec50_ci"]),  # type: ignore[arg-type]
            hill_slope=float(out["hill_slope"]),
            diagnostics=diagnostics,
        )
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Bayesian dose-response failed: {e}")


class RecommendNextCompoundsRequest(BaseModel):
    tested_compounds: List[Dict[str, Any]]
    candidate_pool: List[Dict[str, Any]]
    batch_size: int = 10


@router.post("/analysis/recommend-next-compounds", response_model=List[Dict[str, Any]])
async def recommend_next_compounds_endpoint(
    request: RecommendNextCompoundsRequest,
) -> List[Dict[str, Any]]:
    """Bayesian optimization for next-best compounds."""
    try:
        return recommend_next_compounds(
            tested_compounds=request.tested_compounds,
            candidate_pool=request.candidate_pool,
            batch_size=request.batch_size,
        )
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Recommend-next-compounds failed: {e}")


__all__ = ["router"]


