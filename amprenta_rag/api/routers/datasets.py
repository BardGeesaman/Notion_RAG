from datetime import datetime
"""
API router for Datasets.
"""

import json
import re
from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, Response
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import AnnotationCreate, Dataset, DatasetCreate, DatasetUpdate
from amprenta_rag.api.services import datasets as dataset_service
from amprenta_rag.database.models import Note
from amprenta_rag.models.domain import OmicsType
from amprenta_rag.notebooks import generate_dataset_notebook
from amprenta_rag.utils.optimistic_lock import ConflictError, update_with_lock
from amprenta_rag.utils.uuid_utils import ensure_uuid

router = APIRouter()


@router.post("/", response_model=Dataset, status_code=201)
async def create_dataset(
    dataset: DatasetCreate,
    db: Session = Depends(get_database_session),
) -> Dataset:
    """Create a new dataset."""
    return cast(Dataset, dataset_service.create_dataset(db, dataset))


@router.get("/", response_model=List[Dataset])
async def list_datasets(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    omics_type: Optional[OmicsType] = Query(None, description="Filter by omics type"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    experiment_id: Optional[UUID] = Query(None, description="Filter by experiment ID"),
    db: Session = Depends(get_database_session),
) -> List[Dataset]:
    """List all datasets."""
    omics_type_str = omics_type.value if omics_type else None
    return cast(
        List[Dataset],
        dataset_service.get_datasets(
            db,
            skip=skip,
            limit=limit,
            name_filter=name,
            omics_type=omics_type_str,
            program_id=program_id,
            experiment_id=experiment_id,
        ),
    )


@router.get("/{dataset_id}", response_model=Dataset)
async def get_dataset(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
) -> Dataset:
    """Get a dataset by ID."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")
    return cast(Dataset, dataset)


@router.post("/{dataset_id}/annotations", summary="Add annotation to dataset")
async def add_dataset_annotation(
    dataset_id: UUID,
    annotation: AnnotationCreate,
    db: Session = Depends(get_database_session),
) -> Dict[str, Any]:
    """Add a note/annotation to a dataset."""
    dataset = dataset_service.get_dataset(db, dataset_id)
    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    note = Note(
        entity_type="dataset",
        entity_id=cast(Any, ensure_uuid(dataset_id)),
        annotation_type=annotation.annotation_type,
        content=annotation.text,
    )
    db.add(note)
    db.commit()
    db.refresh(note)

    created_at_val = getattr(note, "created_at", None)
    return {
        "id": str(note.id),
        "entity_type": note.entity_type,
        "entity_id": str(note.entity_id),
        "text": note.content,
        "annotation_type": note.annotation_type,
        "created_at": created_at_val.isoformat() if isinstance(created_at_val, datetime) else None,
    }


@router.get("/{dataset_id}/notebook")
async def download_dataset_notebook(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
) -> Response:
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
) -> Dataset:
    """
    Update a dataset using optimistic locking.

    Clients must provide the current `version` of the dataset. If the stored
    version differs, the request fails with HTTP 409 to prevent lost updates.
    """
    db_dataset = dataset_service.get_dataset(db, dataset_id)
    if not db_dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    update_data = dataset.model_dump(exclude_unset=True)
    expected_version = cast(int, update_data.pop("version"))
    program_ids = update_data.pop("program_ids", None)
    experiment_ids = update_data.pop("experiment_ids", None)

    # Handle omics_type enum conversion
    if "omics_type" in update_data and update_data["omics_type"]:
        update_data["omics_type"] = update_data["omics_type"].value

    # Update program relationships if provided
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program

        db_dataset.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_dataset.programs.append(program)

    # Update experiment relationships if provided
    if experiment_ids is not None:
        from amprenta_rag.api.services.experiments import get_experiment

        db_dataset.experiments.clear()
        for experiment_id in experiment_ids:
            experiment = get_experiment(db, experiment_id)
            if experiment:
                db_dataset.experiments.append(experiment)

    try:
        updated = update_with_lock(db_dataset, cast(Dict[str, Any], update_data), expected_version, db)
    except ConflictError:
        raise HTTPException(status_code=409, detail="Version conflict")

    db.refresh(updated)
    return cast(Dataset, updated)


@router.delete("/{dataset_id}", status_code=204)
async def delete_dataset(
    dataset_id: UUID,
    db: Session = Depends(get_database_session),
) -> None:
    """Delete a dataset."""
    success = dataset_service.delete_dataset(db, dataset_id)
    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")

