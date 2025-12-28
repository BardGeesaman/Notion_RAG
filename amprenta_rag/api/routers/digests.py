"""Digest schedules API (weekly executive digests)."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user_with_company, get_database_session
from amprenta_rag.automation.digest_scheduler import DigestScheduler
from amprenta_rag.database.models import DigestSchedule, Program, User


router = APIRouter(prefix="/digests", tags=["Digests"])


class DigestScheduleCreateRequest(BaseModel):
    program_id: UUID
    notebook_path: str
    schedule_cron: str = Field(..., description="Cron expression (crontab format)")
    recipients: List[str] = Field(default_factory=list)
    enabled: bool = True


class DigestScheduleUpdateRequest(BaseModel):
    notebook_path: Optional[str] = None
    schedule_cron: Optional[str] = None
    recipients: Optional[List[str]] = None
    enabled: Optional[bool] = None


class DigestScheduleResponse(BaseModel):
    id: UUID
    program_id: UUID
    notebook_path: str
    schedule_cron: str
    recipients: List[str]
    last_run_at: Optional[datetime] = None
    last_status: Optional[str] = None
    enabled: bool

    model_config = ConfigDict(from_attributes=True)


@router.post("", response_model=DigestScheduleResponse, status_code=201)
def create_digest_schedule(
    payload: DigestScheduleCreateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> DigestScheduleResponse:
    _ = current_user
    prog = db.query(Program).filter(Program.id == payload.program_id).first()
    if not prog:
        raise HTTPException(status_code=404, detail="Program not found")

    ds = DigestSchedule(
        program_id=payload.program_id,
        notebook_path=payload.notebook_path.strip(),
        schedule_cron=payload.schedule_cron.strip(),
        recipients=list(payload.recipients or []),
        enabled=bool(payload.enabled),
    )
    db.add(ds)
    db.commit()
    db.refresh(ds)

    if ds.enabled:
        try:
            DigestScheduler().add_digest_schedule(UUID(str(ds.id)))
        except Exception:
            # Don't fail API create if scheduler isn't running in this process
            pass

    return DigestScheduleResponse.model_validate(ds)


@router.get("", response_model=List[DigestScheduleResponse])
def list_digest_schedules(
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> List[DigestScheduleResponse]:
    _ = current_user
    items = db.query(DigestSchedule).order_by(DigestSchedule.enabled.desc()).all()
    return [DigestScheduleResponse.model_validate(x) for x in items]


@router.put("/{schedule_id}", response_model=DigestScheduleResponse)
def update_digest_schedule(
    schedule_id: UUID,
    payload: DigestScheduleUpdateRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> DigestScheduleResponse:
    _ = current_user
    ds = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
    if not ds:
        raise HTTPException(status_code=404, detail="Digest schedule not found")

    for k, v in payload.model_dump(exclude_unset=True).items():
        setattr(ds, k, v)
    db.add(ds)
    db.commit()
    db.refresh(ds)

    try:
        sched = DigestScheduler()
        if ds.enabled:
            sched.add_digest_schedule(schedule_id)
        else:
            sched.remove_digest_schedule(schedule_id)
    except Exception:
        pass

    return DigestScheduleResponse.model_validate(ds)


@router.delete("/{schedule_id}", response_model=dict)
def delete_digest_schedule(
    schedule_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    _ = current_user
    ds = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
    if not ds:
        raise HTTPException(status_code=404, detail="Digest schedule not found")
    db.delete(ds)
    db.commit()

    try:
        DigestScheduler().remove_digest_schedule(schedule_id)
    except Exception:
        pass

    return {"deleted": True}


@router.post("/{schedule_id}/run", response_model=dict)
def run_digest_schedule(
    schedule_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    _ = current_user
    ds = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
    if not ds:
        raise HTTPException(status_code=404, detail="Digest schedule not found")

    res = DigestScheduler().run_digest(schedule_id)
    return {"job": "ran", "status": res.status, "output_html": str(res.output_html) if res.output_html else None}


@router.get("/{schedule_id}/output")
def get_last_digest_output(
    schedule_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> FileResponse:
    _ = current_user
    ds = db.query(DigestSchedule).filter(DigestSchedule.id == schedule_id).first()
    if not ds:
        raise HTTPException(status_code=404, detail="Digest schedule not found")

    root = Path(__file__).resolve().parents[3]  # .../RAG
    html_path = root / "data" / "digests" / str(schedule_id) / "latest.html"
    if not html_path.exists():
        raise HTTPException(status_code=404, detail="No digest output available yet")

    return FileResponse(path=str(html_path), media_type="text/html", filename=f"digest_{schedule_id}.html")


__all__ = ["router"]


