"""Batch effect correction endpoints."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from amprenta_rag.analysis.batch_correction import correct_batch_effects

router = APIRouter()


class BatchCorrectRequest(BaseModel):
    datasets: List[UUID]
    batch_map: Dict[str, Any]
    method: str = "combat"


class BatchCorrectResponse(BaseModel):
    corrected_dataset_id: Optional[UUID] = None
    stats: Dict[str, Any]
    preview: List[Dict[str, Any]]


@router.post("/analysis/batch-correct", response_model=BatchCorrectResponse)
def batch_correct(payload: BatchCorrectRequest) -> BatchCorrectResponse:
    try:
        out = correct_batch_effects(
            datasets=payload.datasets,
            batch_map=payload.batch_map,
            method=payload.method,
        )
        df = out["corrected_df"]
        preview = df.reset_index().head(20).to_dict(orient="records")
        return BatchCorrectResponse(
            corrected_dataset_id=out.get("corrected_dataset_id"),
            stats=out["stats"],
            preview=preview,
        )
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Batch correction failed: {e}")


__all__ = ["router"]


