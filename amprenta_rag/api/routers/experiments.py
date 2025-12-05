"""
API router for Experiments.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.api.schemas import Experiment, ExperimentCreate, ExperimentUpdate
from amprenta_rag.api.services import experiments as experiment_service

router = APIRouter()


@router.post("/", response_model=Experiment, status_code=201)
async def create_experiment(
    experiment: ExperimentCreate,
    db: Session = Depends(get_database_session),
):
    """Create a new experiment."""
    return experiment_service.create_experiment(db, experiment)


@router.get("/", response_model=List[Experiment])
async def list_experiments(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    name: Optional[str] = Query(None, description="Filter by name (partial match)"),
    program_id: Optional[UUID] = Query(None, description="Filter by program ID"),
    db: Session = Depends(get_database_session),
):
    """List all experiments."""
    return experiment_service.get_experiments(
        db, skip=skip, limit=limit, name_filter=name, program_id=program_id
    )


@router.get("/{experiment_id}", response_model=Experiment)
async def get_experiment(
    experiment_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get an experiment by ID."""
    experiment = experiment_service.get_experiment(db, experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")
    return experiment


@router.patch("/{experiment_id}", response_model=Experiment)
async def update_experiment(
    experiment_id: UUID,
    experiment: ExperimentUpdate,
    db: Session = Depends(get_database_session),
):
    """Update an experiment."""
    updated = experiment_service.update_experiment(db, experiment_id, experiment)
    if not updated:
        raise HTTPException(status_code=404, detail="Experiment not found")
    return updated


@router.delete("/{experiment_id}", status_code=204)
async def delete_experiment(
    experiment_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Delete an experiment."""
    success = experiment_service.delete_experiment(db, experiment_id)
    if not success:
        raise HTTPException(status_code=404, detail="Experiment not found")

