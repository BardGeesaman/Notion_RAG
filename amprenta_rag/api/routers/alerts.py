from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException

from amprenta_rag.api.dependencies import get_optional_user_id
from amprenta_rag.api.schemas import (
    AlertBatchRequest,
    AlertBatchResponse,
    AlertCheckRequest,
    AlertCheckResponse,
    FilterInfo,
)
from amprenta_rag.api.services import alert_service
from amprenta_rag.chemistry.alert_checker import StructuralAlertChecker, list_available_filters

router = APIRouter()

# Structural alerts router (new): /api/alerts/*
structural_router = APIRouter(prefix="/alerts", tags=["alerts"])


@router.get("/alerts", response_model=List[dict])
def list_alerts(
    unread_only: bool = True,
    user_id: Optional[UUID] = Depends(get_optional_user_id),
):
    return alert_service.list_alerts(user_id, unread_only)


@router.post("/alerts/{alert_id}/read", response_model=dict)
def mark_alert_read(alert_id: UUID):
    ok = alert_service.mark_read(alert_id)
    if not ok:
        raise HTTPException(status_code=404, detail="Alert not found")
    return {"status": "ok"}


@router.post("/alerts/read-all", response_model=dict)
def mark_all_alerts_read(user_id: Optional[UUID] = Depends(get_optional_user_id)):
    count = alert_service.mark_all_read(user_id)
    return {"updated": count}


@structural_router.post("/check", response_model=AlertCheckResponse)
def check_compound(request: AlertCheckRequest) -> AlertCheckResponse:
    chk = StructuralAlertChecker(filters=request.filters)
    out = chk.check_compound(request.smiles)
    return AlertCheckResponse.model_validate(out)


@structural_router.post("/batch", response_model=AlertBatchResponse)
def check_batch(request: AlertBatchRequest) -> AlertBatchResponse:
    if len(request.smiles_list) > 100:
        raise HTTPException(status_code=400, detail="Max 100 SMILES per request")
    chk = StructuralAlertChecker(filters=request.filters)
    results = chk.check_batch(request.smiles_list, n_jobs=4)
    clean = sum(1 for r in results if isinstance(r, dict) and r.get("is_clean") is True)
    flagged = sum(1 for r in results if isinstance(r, dict) and r.get("is_clean") is False)
    return AlertBatchResponse(
        results=[AlertCheckResponse.model_validate(r) for r in results],
        total_checked=len(results),
        clean_count=int(clean),
        flagged_count=int(flagged),
    )


@structural_router.get("/filters", response_model=List[FilterInfo])
def list_filters() -> List[FilterInfo]:
    infos = list_available_filters()
    return [FilterInfo.model_validate(x) for x in infos]


__all__ = ["router", "structural_router"]


