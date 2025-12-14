"""
API router for Datasets.
"""

import json
import re
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, Response
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import AnnotationCreate, Dataset, DatasetCreate, DatasetUpdate
from amprenta_rag.api.services import datasets as dataset_service
from amprenta_rag.database.models import Note
from amprenta_rag.models.domain import OmicsType
from amprenta_rag.notebooks import generate_dataset_notebook

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


@router.post("/{dataset_id}/annotations", summary="Add annotation to dataset")
async def add_dataset_annotation(
    dataset_id: UUID,
    annotation: AnnotationCreate,
    db: Session = Depends(get_database_session),
):
    """Add a note/annotation to a dataset."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    note = Note(
        entity_type="dataset",
        entity_id=dataset_id,
        annotation_type=annotation.annotation_type,
        content=annotation.text,
    )
    db.add(note)
    db.commit()
    db.refresh(note)

    return {
        "id": str(note.id),
        "entity_type": note.entity_type,
        "entity_id": str(note.entity_id),
        "text": note.content,
        "annotation_type": note.annotation_type,
        "created_at": note.created_at.isoformat() if getattr(note, "created_at", None) else None,
    }


@router.get("/{dataset_id}/notebook")
async def download_dataset_notebook(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
):
    """
    Download an nbformat v4 notebook (JSON) for exploring a dataset.
    """
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    nb = generate_dataset_notebook(str(dataset_id))

    dataset_name = getattr(dataset, "name", None) or f"dataset_{str(dataset_id)[:8]}"
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", dataset_name).strip("_") or f"dataset_{str(dataset_id)[:8]}"

    headers = {"Content-Disposition": f'attachment; filename="{safe_name}.ipynb"'}
    return Response(content=json.dumps(nb), media_type="application/json", headers=headers)


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

