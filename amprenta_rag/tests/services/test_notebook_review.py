from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any
from uuid import UUID, uuid4

from amprenta_rag.services.notebook_review import NotebookReviewService, _compute_signature


class _Query:
    def __init__(self, first_obj: Any = None):
        self._first = first_obj

    def filter(self, *_args: Any, **_kwargs: Any) -> "_Query":
        return self

    def order_by(self, *_args: Any, **_kwargs: Any) -> "_Query":
        return self

    def first(self) -> Any:
        return self._first


class _DB:
    def __init__(self):
        self.obj = None
        self.commits = 0

    def add(self, obj: Any) -> None:
        self.obj = obj

    def commit(self) -> None:
        self.commits += 1

    def refresh(self, _obj: Any) -> None:
        return None

    def query(self, _model: Any) -> _Query:
        return _Query(first_obj=self.obj)


@dataclass
class _Review:
    id: UUID
    notebook_path: str
    version_hash: str
    reviewer_id: UUID | None = None
    status: str = "pending"
    comments: str | None = None
    signature: str | None = None
    reviewed_at: datetime | None = None
    created_at: datetime | None = None


def test_request_review_creates_pending(monkeypatch):
    os.environ["NOTEBOOK_REVIEW_SECRET"] = "testsecret"
    db = _DB()

    # Avoid filesystem/registry access
    monkeypatch.setattr("amprenta_rag.services.notebook_review._resolve_and_hash_notebook", lambda p: "h" * 64)

    # Use a lightweight stand-in model object that matches attribute usage
    monkeypatch.setattr(
        "amprenta_rag.services.notebook_review.NotebookReview",
        lambda **kwargs: _Review(id=uuid4(), **kwargs),
    )

    svc = NotebookReviewService(db)
    rid = svc.request_review("deploy/jupyterhub/templates/dose_response.ipynb", requester_id=uuid4())
    assert isinstance(rid, UUID)
    assert db.obj.status == "pending"
    assert db.obj.version_hash == "h" * 64


def test_submit_approval_generates_signature(monkeypatch):
    os.environ["NOTEBOOK_REVIEW_SECRET"] = "testsecret"
    db = _DB()
    review_id = uuid4()
    db.obj = _Review(id=review_id, notebook_path="x.ipynb", version_hash="a" * 64)

    svc = NotebookReviewService(db)
    sig = svc.submit_review(review_id, reviewer_id=uuid4(), status="approved", comments="ok")
    assert isinstance(sig, str) and len(sig) == 64
    assert db.obj.status == "approved"
    assert db.obj.signature == sig
    assert db.obj.reviewed_at is not None


def test_verify_signature_valid():
    os.environ["NOTEBOOK_REVIEW_SECRET"] = "testsecret"
    db = _DB()
    rid = uuid4()
    reviewer = uuid4()
    reviewed_at = datetime.now(timezone.utc)
    version_hash = "b" * 64
    sig = _compute_signature(reviewer, version_hash, reviewed_at)
    db.obj = _Review(
        id=rid,
        notebook_path="x.ipynb",
        version_hash=version_hash,
        reviewer_id=reviewer,
        status="approved",
        signature=sig,
        reviewed_at=reviewed_at,
        created_at=reviewed_at,
    )

    svc = NotebookReviewService(db)
    assert svc.verify_signature(rid) is True


def test_verify_signature_tampered():
    os.environ["NOTEBOOK_REVIEW_SECRET"] = "testsecret"
    db = _DB()
    rid = uuid4()
    reviewer = uuid4()
    reviewed_at = datetime.now(timezone.utc)
    version_hash = "c" * 64
    sig = _compute_signature(reviewer, version_hash, reviewed_at)
    db.obj = _Review(
        id=rid,
        notebook_path="x.ipynb",
        version_hash=version_hash,
        reviewer_id=reviewer,
        status="approved",
        signature=sig + "00",  # tamper
        reviewed_at=reviewed_at,
        created_at=reviewed_at,
    )

    svc = NotebookReviewService(db)
    assert svc.verify_signature(rid) is False


