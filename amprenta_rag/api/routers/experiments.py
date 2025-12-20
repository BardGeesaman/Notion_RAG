 # mypy: ignore-errors
"""
API router for Experiments.
"""

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import AnnotationCreate, Experiment, ExperimentCreate, ExperimentUpdate
from amprenta_rag.api.services import experiments as experiment_service
from amprenta_rag.database.models import Note

router = APIRouter()


@router.post("/", response_model=Experiment, status_code=201)
async def create_experiment(
    experiment: ExperimentCreate,
    db: Session = Depends(get_database_session),
) -> Experiment:
    """Create a new experiment."""
    return experiment_service.create_experiment(db, experiment)


@router.get("/", response_model=List[Experiment])
async def list_experiments(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    db: Session = Depends(get_database_session),
) -> List[Experiment]:
    """List all experiments."""
    return experiment_service.get_experiments(
        db, skip=skip, limit=limit, name_filter=name, program_id=program_id
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
    return experiment


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
        entity_id=experiment_id,  # type: ignore[arg-type]
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


@router.patch("/{experiment_id}", response_model=Experiment)
async def update_experiment(
    experiment_id: UUID,
    experiment: ExperimentUpdate,
    db: Session = Depends(get_database_session),
) -> Experiment:
    """Update an experiment."""
    updated = experiment_service.update_experiment(db, experiment_id, experiment)
    if not updated:
        raise HTTPException(status_code=404, detail="Experiment not found")
    return updated


@router.delete("/{experiment_id}", status_code=204)
async def delete_experiment(
    experiment_id: UUID,
    db: Session = Depends(get_database_session),
) -> None:
    """Delete an experiment."""
    success = experiment_service.delete_experiment(db, experiment_id)
    if not success:
        raise HTTPException(status_code=404, detail="Experiment not found")

