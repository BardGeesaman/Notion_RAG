from __future__ import annotations

from datetime import datetime, timezone, timedelta
from typing import Dict, List, Optional

from sqlalchemy import func, or_
from sqlalchemy.orm import selectinload

from amprenta_rag.database.models import Dataset
from amprenta_rag.database.session import db_session

# pylint: disable=not-callable


def _normalize_dt(dt: Optional[datetime]) -> Optional[datetime]:
    """Ensure datetime is timezone-aware for comparisons."""
    if dt is None:
        return None
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt


def _health_status(last_sync: Optional[datetime]) -> str:
    """Map last sync time to health buckets."""
    last_sync = _normalize_dt(last_sync)
    if not last_sync:
        return "stale"
    now = datetime.now(timezone.utc)
    delta: timedelta = now - last_sync
    if delta.days < 7:
        return "healthy"
    if delta.days > 30:
        return "stale"
    return "warning"


def get_repository_summary() -> List[Dict[str, object]]:
    """Aggregate dataset counts and recency by source (data_origin)."""
    with db_session() as db:
        rows = (
            db.query(
                Dataset.data_origin,
                func.count(Dataset.id).label("count"),
                func.max(Dataset.updated_at).label("last_sync"),
            )
            .group_by(Dataset.data_origin)
            .all()
        )

    summaries: List[Dict[str, object]] = []
    for source, count, last_sync in rows:
        summaries.append(
            {
                "name": source or "Unknown",
                "dataset_count": int(count or 0),
                "last_sync_date": _normalize_dt(last_sync),
                "health_status": _health_status(last_sync),
            }
        )
    return summaries


def get_catalog_datasets(
    source: Optional[str] = None,
    search: Optional[str] = None,
    limit: int = 50,
    offset: int = 0,
) -> tuple[List[dict], int]:
    """List datasets with optional source/search filter and pagination."""
    with db_session() as db:
        query = db.query(Dataset).options(selectinload(Dataset.features))
        if source:
            query = query.filter(Dataset.data_origin == source)
        if search:
            term = f"%{search}%"
            # Dataset model lacks explicit title/accession columns; use name and external_ids text
            query = query.filter(
                or_(
                    Dataset.name.ilike(term),
                    Dataset.description.ilike(term),
                )
            )

        total = query.count()
        results: List[Dataset] = (
            query.order_by(Dataset.created_at.desc()).offset(offset).limit(limit).all()
        )

        datasets: List[dict] = []
        for d in results:
            ext = getattr(d, "external_ids", {}) or {}
            accession = ext.get("accession") or ext.get("study_id") or ext.get("id")
            datasets.append(
                {
                    "id": d.id,
                    "accession": accession,
                    "title": getattr(d, "name", None),
                    "source": getattr(d, "data_origin", None),
                    "created_at": d.created_at,
                    "status": getattr(d, "ingestion_status", None),
                    "feature_count": len(getattr(d, "features", []) or []),
                }
            )

    return datasets, total

