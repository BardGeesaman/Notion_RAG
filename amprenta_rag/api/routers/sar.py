from __future__ import annotations

from typing import List

from fastapi import APIRouter, Query

from amprenta_rag.api import schemas
from amprenta_rag.api.services import sar_data as service

router = APIRouter()


@router.get(
    "/targets",
    summary="List SAR targets",
    response_model=List[schemas.TargetResponse],
)
def list_targets(limit: int = Query(200, ge=1, le=2000)):
    return service.list_targets(limit=limit)


@router.get(
    "/targets/{target}/compounds",
    summary="List compounds (with activity) for a target",
    response_model=List[schemas.CompoundActivityResponse],
)
def get_compounds_by_target(target: str, limit: int = Query(2000, ge=1, le=20000)):
    return service.get_compounds_by_target(target=target, limit=limit)


@router.get(
    "/targets/{target}/cliffs",
    summary="Detect activity cliffs for a target",
    response_model=List[schemas.ActivityCliffResponse],
)
def get_activity_cliffs(
    target: str,
    similarity_threshold: float = Query(0.6, ge=0.0, le=1.0),
    fold_change: float = Query(10.0, ge=1.0, le=1e6),
    limit: int = Query(50, ge=1, le=500),
):
    return service.get_activity_cliffs_for_target(
        target=target,
        similarity_threshold=similarity_threshold,
        fold_change=fold_change,
        limit=limit,
    )


