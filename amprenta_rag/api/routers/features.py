 # mypy: ignore-errors
"""
API router for Features.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import Feature, FeatureCreate, FeatureUpdate
from amprenta_rag.api.services import features as feature_service
from amprenta_rag.models.domain import FeatureType

router = APIRouter()


@router.post("/", response_model=Feature, status_code=201)
async def create_feature(
    feature: FeatureCreate,
    db: Session = Depends(get_database_session),
):
    """Create a new feature."""
    return feature_service.create_feature(db, feature)


@router.get("/", response_model=List[Feature])
async def list_features(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    feature_type: Optional[FeatureType] = Query(None, description="Filter by feature type"),
    dataset_id: Optional[UUID] = Query(None, description="Filter by dataset ID"),
    db: Session = Depends(get_database_session),
):
    """List all features."""
    return feature_service.get_features(
        db,
        skip=skip,
        limit=limit,
        name_filter=name,
        feature_type=feature_type,
        dataset_id=dataset_id,
    )


@router.get("/{feature_id}", response_model=Feature)
async def get_feature(
    feature_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get a feature by ID."""
    feature = feature_service.get_feature(db, feature_id)
    if not feature:
        raise HTTPException(status_code=404, detail="Feature not found")
    return feature


@router.patch("/{feature_id}", response_model=Feature)
async def update_feature(
    feature_id: UUID,
    feature: FeatureUpdate,
    db: Session = Depends(get_database_session),
):
    """Update a feature."""
    updated = feature_service.update_feature(db, feature_id, feature)
    if not updated:
        raise HTTPException(status_code=404, detail="Feature not found")
    return updated


@router.delete("/{feature_id}", status_code=204)
async def delete_feature(
    feature_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Delete a feature."""
    success = feature_service.delete_feature(db, feature_id)
    if not success:
        raise HTTPException(status_code=404, detail="Feature not found")

