"""Sphingolipid analysis endpoints."""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from amprenta_rag.analysis.sphingolipid.scoring import compute_pathway_imbalance

router = APIRouter()


class SphingolipidImbalanceRequest(BaseModel):
    dataset_ids: List[UUID]
    pathway: str = "ceramide"


class SphingolipidImbalanceResponse(BaseModel):
    pathway: str
    method: str
    score: float
    matched: Dict[str, Any]
    stats: Dict[str, Any]


@router.post("/analysis/sphingolipid/imbalance", response_model=SphingolipidImbalanceResponse)
def sphingolipid_imbalance(payload: SphingolipidImbalanceRequest) -> SphingolipidImbalanceResponse:
    try:
        out = compute_pathway_imbalance(payload.dataset_ids, pathway=payload.pathway)
        return SphingolipidImbalanceResponse(**out)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Sphingolipid imbalance failed: {e}")


__all__ = ["router"]


