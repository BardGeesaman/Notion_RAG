"""Notebook review workflow service (signed review cards).

Signature scheme:
  HMAC-SHA256(reviewer_id + version_hash + timestamp + secret)

Timestamp is the persisted reviewed_at in ISO8601 format.
"""

from __future__ import annotations

import hashlib
import hmac
import os
from datetime import datetime, timezone
from typing import Any, Dict
from uuid import UUID

from amprenta_rag.database.models import NotebookReview, ActivityEventType
from amprenta_rag.notebooks.registry import load_registry, resolve_notebook_path
from amprenta_rag.services.activity import log_activity
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Import notebook_diff service for snapshot capture
try:
    from amprenta_rag.services.notebook_diff import capture_snapshot
except ImportError:
    # Handle case where notebook_diff service is not available
    capture_snapshot = None


VALID_STATUSES = {"pending", "approved", "rejected", "changes_requested"}


def _get_secret() -> bytes:
    secret = (os.environ.get("NOTEBOOK_REVIEW_SECRET") or os.environ.get("REVIEW_HMAC_SECRET") or "").strip()
    if not secret:
        raise RuntimeError("Missing NOTEBOOK_REVIEW_SECRET (or REVIEW_HMAC_SECRET)")
    return secret.encode("utf-8")


def _resolve_and_hash_notebook(notebook_path: str) -> str:
    nb_rel = (notebook_path or "").strip()
    if not nb_rel:
        raise ValueError("notebook_path is required")

    for tpl in load_registry():
        if not isinstance(tpl, dict):
            continue
        if str(tpl.get("notebook_path") or "").strip() != nb_rel:
            continue
        p = resolve_notebook_path(tpl)
        data = p.read_bytes()
        return hashlib.sha256(data).hexdigest()
    raise FileNotFoundError("notebook_path not found in registry")


def _compute_signature(reviewer_id: UUID, version_hash: str, reviewed_at: datetime) -> str:
    ts = reviewed_at.astimezone(timezone.utc).isoformat()
    msg = f"{reviewer_id}|{version_hash}|{ts}".encode("utf-8")
    sig = hmac.new(_get_secret(), msg, digestmod=hashlib.sha256).hexdigest()
    return sig


class NotebookReviewService:
    """Service for requesting and submitting notebook reviews."""

    def __init__(self, db_session):
        self.db = db_session

    def request_review(self, notebook_path: str, requester_id: UUID) -> UUID:
        """Create a pending review for the current notebook content."""
        _ = requester_id  # schema does not store requester_id; audit can be handled elsewhere.

        version_hash = _resolve_and_hash_notebook(notebook_path)
        r = NotebookReview(
            notebook_path=notebook_path.strip(),
            version_hash=version_hash,
            status="pending",
            comments=None,
            reviewer_id=None,
            signature=None,
            reviewed_at=None,
            created_at=datetime.now(timezone.utc),
        )
        self.db.add(r)
        self.db.commit()
        self.db.refresh(r)
        
        # Capture notebook snapshot for diff computation
        if capture_snapshot:
            try:
                capture_snapshot(self.db, r.id, notebook_path)
                logger.info(f"Captured snapshot for review {r.id}")
            except Exception as e:
                # Don't fail the review request if snapshot capture fails
                logger.warning(f"Failed to capture snapshot for review {r.id}: {e}")
        else:
            logger.debug("Snapshot capture not available")
        
        return UUID(str(r.id))

    def submit_review(self, review_id: UUID, reviewer_id: UUID, status: str, comments: str = "") -> str:
        """Submit a review decision and return the generated signature."""
        st = (status or "").strip().lower()
        if st not in VALID_STATUSES or st == "pending":
            raise ValueError("status must be approved|rejected|changes_requested")

        r = self.db.query(NotebookReview).filter(NotebookReview.id == review_id).first()
        if not r:
            raise KeyError("review not found")

        reviewed_at = datetime.now(timezone.utc)
        sig = _compute_signature(reviewer_id=reviewer_id, version_hash=str(r.version_hash), reviewed_at=reviewed_at)

        r.reviewer_id = reviewer_id
        r.status = st
        r.comments = comments
        r.reviewed_at = reviewed_at
        r.signature = sig
        self.db.add(r)
        self.db.commit()
        self.db.refresh(r)
        
        # Log activity
        try:
            log_activity(
                event_type=ActivityEventType.NOTEBOOK_REVIEWED,
                target_type="notebook",
                target_id=r.id,
                target_name=r.notebook_path.split('/')[-1] if r.notebook_path else str(r.id),
                actor_id=reviewer_id,
                program_id=None,  # No program context in notebook reviews
                metadata={
                    "verdict": r.status,
                    "notebook_path": r.notebook_path,
                    "version_hash": r.version_hash,
                }
            )
        except Exception as e:
            logger.error(f"Failed to log notebook review activity: {e}")
        
        return sig

    def get_review_status(self, notebook_path: str) -> Dict[str, Any]:
        """Return the latest review status for a notebook."""
        nb = (notebook_path or "").strip()
        if not nb:
            raise ValueError("notebook_path is required")
        r = (
            self.db.query(NotebookReview)
            .filter(NotebookReview.notebook_path == nb)
            .order_by(NotebookReview.created_at.desc())
            .first()
        )
        if not r:
            return {"status": None, "reviewer_id": None, "signature": None, "reviewed_at": None, "version_hash": None}
        return {
            "status": r.status,
            "reviewer_id": str(r.reviewer_id) if r.reviewer_id else None,
            "signature": r.signature,
            "reviewed_at": r.reviewed_at.isoformat() if r.reviewed_at else None,
            "version_hash": r.version_hash,
        }

    def verify_signature(self, review_id: UUID) -> bool:
        """Verify signature integrity for a review."""
        r = self.db.query(NotebookReview).filter(NotebookReview.id == review_id).first()
        if not r:
            raise KeyError("review not found")
        if not r.reviewer_id or not r.reviewed_at or not r.signature:
            return False
        expected = _compute_signature(
            reviewer_id=UUID(str(r.reviewer_id)),
            version_hash=str(r.version_hash),
            reviewed_at=r.reviewed_at,
        )
        return hmac.compare_digest(str(r.signature), expected)


__all__ = ["NotebookReviewService", "VALID_STATUSES"]


