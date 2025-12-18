from __future__ import annotations

from typing import Optional, List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException

from amprenta_rag.api.dependencies import get_optional_user_id
from amprenta_rag.api.services import alert_service

router = APIRouter()


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

