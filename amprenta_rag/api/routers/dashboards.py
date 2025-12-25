"""Pinned dashboards API endpoints (Voila notebook dashboards)."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user_with_company, get_database_session
from amprenta_rag.database.models import PinnedDashboard, Program, User
from amprenta_rag.notebooks.registry import load_registry, resolve_notebook_path


router = APIRouter(prefix="/dashboards", tags=["Dashboards"])


def _validate_notebook_path(nb_path: str) -> None:
    """Ensure notebook_path exists in registry and resolves to an on-disk .ipynb in templates."""
    p = (nb_path or "").strip()
    if not p:
        raise HTTPException(status_code=400, detail="notebook_path is required")
    for tpl in load_registry():
        if not isinstance(tpl, dict):
            continue
        if str(tpl.get("notebook_path") or "").strip() != p:
            continue
        try:
            abs_path = resolve_notebook_path(tpl)
        except Exception:
            raise HTTPException(status_code=400, detail="invalid notebook_path")
        if not abs_path.exists():
            raise HTTPException(status_code=404, detail="notebook not found on disk")
        if abs_path.suffix != ".ipynb":
            raise HTTPException(status_code=400, detail="invalid notebook type")
        return
    raise HTTPException(status_code=404, detail="notebook_path not found in registry")


class PinDashboardRequest(BaseModel):
    notebook_path: str = Field(..., description="Relative notebook path under templates/")
    program_id: UUID
    display_name: Optional[str] = None
    config: Optional[Dict[str, Any]] = None


class PinnedDashboardResponse(BaseModel):
    id: UUID
    notebook_path: str
    program_id: UUID
    pinned_by_id: Optional[UUID] = None
    pinned_at: datetime
    display_name: Optional[str] = None
    config: Optional[Dict[str, Any]] = None

    class Config:
        from_attributes = True


@router.post("/pin", response_model=PinnedDashboardResponse, status_code=201)
def pin_dashboard(
    payload: PinDashboardRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> PinnedDashboardResponse:
    _validate_notebook_path(payload.notebook_path)

    prog = db.query(Program).filter(Program.id == payload.program_id).first()
    if not prog:
        raise HTTPException(status_code=404, detail="Program not found")

    pin = PinnedDashboard(
        notebook_path=payload.notebook_path.strip(),
        program_id=payload.program_id,
        pinned_by_id=getattr(current_user, "id", None),
        pinned_at=datetime.now(timezone.utc),
        display_name=(payload.display_name.strip() if payload.display_name else None),
        config=payload.config,
    )
    db.add(pin)
    db.commit()
    db.refresh(pin)
    return PinnedDashboardResponse.model_validate(pin)


@router.delete("/pin/{pin_id}", response_model=dict)
def unpin_dashboard(
    pin_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    pin = db.query(PinnedDashboard).filter(PinnedDashboard.id == pin_id).first()
    if not pin:
        raise HTTPException(status_code=404, detail="Pinned dashboard not found")
    # Only owner or admin can remove
    if pin.pinned_by_id and pin.pinned_by_id != getattr(current_user, "id", None):
        if (getattr(current_user, "role", None) or "").lower() != "admin":
            raise HTTPException(status_code=403, detail="Forbidden")
    db.delete(pin)
    db.commit()
    return {"unpinned": True}


@router.get("/program/{program_id}", response_model=List[PinnedDashboardResponse])
def list_program_dashboards(
    program_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> List[PinnedDashboardResponse]:
    _ = current_user  # auth enforced
    pins = (
        db.query(PinnedDashboard)
        .filter(PinnedDashboard.program_id == program_id)
        .order_by(PinnedDashboard.pinned_at.desc())
        .all()
    )
    return [PinnedDashboardResponse.model_validate(p) for p in pins]


@router.get("/mine", response_model=List[PinnedDashboardResponse])
def list_my_dashboards(
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> List[PinnedDashboardResponse]:
    uid = getattr(current_user, "id", None)
    pins = (
        db.query(PinnedDashboard)
        .filter(PinnedDashboard.pinned_by_id == uid)
        .order_by(PinnedDashboard.pinned_at.desc())
        .all()
    )
    return [PinnedDashboardResponse.model_validate(p) for p in pins]


__all__ = ["router"]



