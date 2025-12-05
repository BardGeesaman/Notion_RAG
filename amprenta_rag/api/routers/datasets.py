"""
API router for Datasets.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import Dataset, DatasetCreate, DatasetUpdate
from amprenta_rag.api.services import datasets as dataset_service
from amprenta_rag.models.domain import OmicsType

router = APIRouter()


@router.post("/", response_model=Dataset, status_code=201)
async def create_dataset(
    dataset: DatasetCreate,
    db: Session = Depends(get_database_session),
):
    """Create a new dataset."""
    return dataset_service.create_dataset(db, dataset)


@router.get("/", response_model=List[Dataset])
async def list_datasets(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    omics_type: Optional[OmicsType] = Query(None, description="Filter by omics type"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    experiment_id: Optional[UUID] = Query(None, description="Filter by experiment ID"),
    db: Session = Depends(get_database_session),
):
    """List all datasets."""
    omics_type_str = omics_type.value if omics_type else None
    return dataset_service.get_datasets(
        db,
        skip=skip,
        limit=limit,
        name_filter=name,
        omics_type=omics_type_str,
        program_id=program_id,
        experiment_id=experiment_id,
    )


@router.get("/{dataset_id}", response_model=Dataset)
async def get_dataset(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get a dataset by ID."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return dataset


@router.patch("/{dataset_id}", response_model=Dataset)
async def update_dataset(
    dataset_id: UUID,
    dataset: DatasetUpdate,
    db: Session = Depends(get_database_session),
):
    """Update a dataset."""
    updated = dataset_service.update_dataset(db, dataset_id, dataset)
    if not updated:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return updated


@router.delete("/{dataset_id}", status_code=204)
async def delete_dataset(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Delete a dataset."""
    success = dataset_service.delete_dataset(db, dataset_id)
    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")

