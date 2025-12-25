"""Compound-target network API endpoints (GraphEdge-backed)."""

from __future__ import annotations

from typing import Dict, List
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.schemas import (
    CompoundTargetNetworkRequest,
    CompoundTargetNetworkResponse,
)
from amprenta_rag.graph.compound_target_network import CompoundTargetNetworkService


router = APIRouter(prefix="/compound-target", tags=["compound-target-network"])


@router.post("", response_model=CompoundTargetNetworkResponse)
def build_network(payload: CompoundTargetNetworkRequest) -> CompoundTargetNetworkResponse:
    svc = CompoundTargetNetworkService()
    try:
        out = svc.get_compound_target_network(
            compound_ids=payload.compound_ids,
            target_ids=payload.target_ids,
            filters=payload.filters,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=str(e))
    return CompoundTargetNetworkResponse.model_validate(out)


@router.get("/expand/compound/{compound_id}")
def expand_from_compound(compound_id: UUID) -> Dict[str, List[str]]:
    svc = CompoundTargetNetworkService()
    try:
        ids = svc.expand_from_compound(compound_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"target_ids": [str(x) for x in ids]}


@router.get("/expand/target/{target_id}")
def expand_from_target(target_id: UUID) -> Dict[str, List[str]]:
    svc = CompoundTargetNetworkService()
    try:
        ids = svc.expand_from_target(target_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"compound_ids": [str(x) for x in ids]}


__all__ = ["router"]


