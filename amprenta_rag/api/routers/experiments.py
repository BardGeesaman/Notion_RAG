from datetime import datetime
"""
API router for Experiments.
"""

from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import AnnotationCreate, Experiment, ExperimentCreate, ExperimentUpdate
from amprenta_rag.api.services import experiments as experiment_service
from amprenta_rag.database.models import Note
from amprenta_rag.utils.optimistic_lock import ConflictError, update_with_lock
from amprenta_rag.utils.uuid_utils import ensure_uuid

router = APIRouter()


@router.post("/", response_model=Experiment, status_code=201)
async def create_experiment(
    experiment: ExperimentCreate,
    db: Session = Depends(get_database_session),
) -> Experiment:
    """Create a new experiment."""
    return cast(Experiment, experiment_service.create_experiment(db, experiment))


@router.get("/", response_model=List[Experiment])
async def list_experiments(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    db: Session = Depends(get_database_session),
) -> List[Experiment]:
    """List all experiments."""
    return cast(
        List[Experiment],
        experiment_service.get_experiments(
            db, skip=skip, limit=limit, name_filter=name, program_id=program_id
        ),
    )


@router.get("/{experiment_id}", response_model=Experiment)
async def get_experiment(
    experiment_id: UUID,
    db: Session = Depends(get_database_session),
) -> Experiment:
    """Get an experiment by ID."""
    experiment = experiment_service.get_experiment(db, experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")
    return cast(Experiment, experiment)


@router.post("/{experiment_id}/annotations", summary="Add annotation to experiment")
async def add_experiment_annotation(
    experiment_id: UUID,
    annotation: AnnotationCreate,
    db: Session = Depends(get_database_session),
) -> Dict[str, Any]:
    """Add a note/annotation to an experiment."""
    experiment = experiment_service.get_experiment(db, experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    note = Note(
        entity_type="experiment",
        entity_id=cast(Any, ensure_uuid(experiment_id)),
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


@router.patch("/{experiment_id}", response_model=Experiment)
async def update_experiment(
    experiment_id: UUID,
    experiment: ExperimentUpdate,
    db: Session = Depends(get_database_session),
) -> Experiment:
    """
    Update an experiment using optimistic locking.

    Clients must provide the current `version` of the experiment. If the stored
    version differs, the request fails with HTTP 409 to prevent lost updates.
    """
    db_experiment = experiment_service.get_experiment(db, experiment_id)
    if not db_experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    update_data = experiment.model_dump(exclude_unset=True)
    expected_version = cast(int, update_data.pop("version"))
    program_ids = update_data.pop("program_ids", None)

    # Update program relationships (if provided) in the same transaction/version bump.
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program

        db_experiment.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_experiment.programs.append(program)

    try:
        updated = update_with_lock(db_experiment, cast(Dict[str, Any], update_data), expected_version, db)
    except ConflictError:
        raise HTTPException(status_code=409, detail="Version conflict")

    db.refresh(updated)
    return cast(Experiment, updated)


@router.delete("/{experiment_id}", status_code=204)
async def delete_experiment(
    experiment_id: UUID,
    db: Session = Depends(get_database_session),
) -> None:
    """Delete an experiment."""
    success = experiment_service.delete_experiment(db, experiment_id)
    if not success:
        raise HTTPException(status_code=404, detail="Experiment not found")

