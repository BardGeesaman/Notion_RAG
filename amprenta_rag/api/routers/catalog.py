from __future__ import annotations

from fastapi import APIRouter, Query

from amprenta_rag.api.schemas import RepositorySummary, CatalogDataset
from amprenta_rag.api.services import catalog_service

router = APIRouter()


@router.get("/catalog/summary", response_model=dict)
def get_catalog_summary():
    """Return repository-level counts and recency."""
    summaries = catalog_service.get_repository_summary()
    return {"repositories": [RepositorySummary(**s) for s in summaries]}


@router.get("/catalog/datasets", response_model=dict)
def list_catalog_datasets(
    source: str | None = Query(default=None),
    search: str | None = Query(default=None),
    limit: int = Query(default=50, ge=1, le=200),
    offset: int = Query(default=0, ge=0),
):
    """List catalog datasets with optional source/search filters and pagination."""
    datasets_raw, total = catalog_service.get_catalog_datasets(source, search, limit, offset)
    datasets = [CatalogDataset(**d) for d in datasets_raw]
    return {"datasets": datasets, "total": total}

