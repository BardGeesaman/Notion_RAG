"""Notebook review workflow API."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user_with_company, get_database_session
from amprenta_rag.database.models import NotebookReview, User
from amprenta_rag.services.notebook_review import NotebookReviewService, VALID_STATUSES


router = APIRouter(prefix="/reviews", tags=["Reviews"])


class ReviewRequest(BaseModel):
    notebook_path: str


class ReviewSubmitRequest(BaseModel):
    status: str = Field(..., description="approved|rejected|changes_requested")
    comments: Optional[str] = ""


class ReviewResponse(BaseModel):
    id: UUID
    notebook_path: str
    version_hash: str
    reviewer_id: Optional[UUID] = None
    status: str
    comments: Optional[str] = None
    signature: Optional[str] = None
    reviewed_at: Optional[datetime] = None
    created_at: datetime

    class Config:
        from_attributes = True


@router.post("", response_model=dict, status_code=201)
def request_review(
    payload: ReviewRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    svc = NotebookReviewService(db)
    rid = svc.request_review(payload.notebook_path, requester_id=UUID(str(current_user.id)))
    return {"review_id": str(rid)}


@router.put("/{review_id}", response_model=dict)
def submit_review(
    review_id: UUID,
    payload: ReviewSubmitRequest,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    status = (payload.status or "").strip().lower()
    if status not in VALID_STATUSES or status == "pending":
        raise HTTPException(status_code=400, detail="Invalid status")
    svc = NotebookReviewService(db)
    try:
        sig = svc.submit_review(review_id, reviewer_id=UUID(str(current_user.id)), status=status, comments=payload.comments or "")
    except KeyError:
        raise HTTPException(status_code=404, detail="Review not found")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"signature": sig}


@router.get("/notebook/{path:path}", response_model=Dict[str, Any])
def get_notebook_review_status(
    path: str,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> Dict[str, Any]:
    _ = current_user
    svc = NotebookReviewService(db)
    return svc.get_review_status(path)


@router.get("/pending", response_model=List[ReviewResponse])
def list_pending_reviews(
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> List[ReviewResponse]:
    _ = current_user
    # Pending queue: include unassigned or assigned-to-me.
    uid = UUID(str(current_user.id))
    items = (
        db.query(NotebookReview)
        .filter(NotebookReview.status == "pending")
        .order_by(NotebookReview.created_at.desc())
        .all()
    )
    out = []
    for r in items:
        if r.reviewer_id is None or UUID(str(r.reviewer_id)) == uid:
            out.append(ReviewResponse.model_validate(r))
    return out


@router.get("/{review_id}/verify", response_model=dict)
def verify_review_signature(
    review_id: UUID,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user_with_company),
) -> dict:
    _ = current_user
    svc = NotebookReviewService(db)
    try:
        ok = svc.verify_signature(review_id)
    except KeyError:
        raise HTTPException(status_code=404, detail="Review not found")
    return {"valid": bool(ok)}


__all__ = ["router"]


