"""
CRUD services for Experiments.
"""

from typing import List, Optional, Any, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Experiment as ExperimentModel, ActivityEventType
from amprenta_rag.api.schemas import ExperimentCreate, ExperimentUpdate
from amprenta_rag.services.activity import log_activity
from amprenta_rag.logging_utils import get_logger
import uuid

logger = get_logger(__name__)


def create_experiment(db: Session, experiment: ExperimentCreate) -> ExperimentModel:
    """Create a new experiment."""
    db_experiment = ExperimentModel(
        id=cast(Any, uuid.uuid4()),
        name=experiment.name,
        type=experiment.type,
        description=experiment.description,
        disease=cast(Any, experiment.disease if experiment.disease is not None else []),
        matrix=cast(Any, experiment.matrix if experiment.matrix is not None else []),
        model_systems=cast(Any, experiment.model_systems if experiment.model_systems is not None else []),
    )

    # Add program relationships
    if experiment.program_ids:
        from amprenta_rag.api.services.programs import get_program
        for program_id in experiment.program_ids:
            program = get_program(db, program_id)
            if program:
                db_experiment.programs.append(program)

    db.add(db_experiment)
    db.commit()
    db.refresh(db_experiment)
    
    # Log activity
    try:
        # Get program_id from first associated program if available
        program_id = db_experiment.programs[0].id if db_experiment.programs else None
        
        log_activity(
            event_type=ActivityEventType.EXPERIMENT_CREATED,
            target_type="experiment",
            target_id=db_experiment.id,
            target_name=db_experiment.name,
            actor_id=None,  # NOTE: requires service layer refactor for user context
            program_id=program_id,
            metadata={
                "experiment_type": db_experiment.type,
                "description": db_experiment.description,
            }
        )
    except Exception as e:
        logger.error(f"Failed to log experiment creation activity: {e}")
    
    return db_experiment


def get_experiment(db: Session, experiment_id: UUID) -> Optional[ExperimentModel]:
    """Get an experiment by ID."""
    return db.query(ExperimentModel).filter(ExperimentModel.id == experiment_id).first()


def get_experiments(
    db: Session,
    skip: int = 0,
    limit: int = 100,
    name_filter: Optional[str] = None,
    program_id: Optional[UUID] = None,
) -> List[ExperimentModel]:
    """Get all experiments with optional filtering."""
    query = db.query(ExperimentModel)

    if name_filter:
        query = query.filter(ExperimentModel.name.ilike(f"%{name_filter}%"))

    if program_id:
        query = query.filter(ExperimentModel.programs.any(id=program_id))

    return query.offset(skip).limit(limit).all()


def update_experiment(
    db: Session,
    experiment_id: UUID,
    experiment: ExperimentUpdate,
) -> Optional[ExperimentModel]:
    """Update an experiment."""
    db_experiment = get_experiment(db, experiment_id)
    if not db_experiment:
        return None

    update_data = experiment.model_dump(exclude_unset=True)
    program_ids = update_data.pop("program_ids", None)

    for field, value in update_data.items():
        setattr(db_experiment, field, value)

    # Update program relationships if provided
    if program_ids is not None:
        from amprenta_rag.api.services.programs import get_program
        db_experiment.programs.clear()
        for program_id in program_ids:
            program = get_program(db, program_id)
            if program:
                db_experiment.programs.append(program)

    db.commit()
    db.refresh(db_experiment)
    return db_experiment


def delete_experiment(db: Session, experiment_id: UUID) -> bool:
    """Delete an experiment."""
    db_experiment = get_experiment(db, experiment_id)
    if not db_experiment:
        return False

    db.delete(db_experiment)
    db.commit()
    return True

