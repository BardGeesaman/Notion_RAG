"""
Experiment CRUD operations.

Split out of `amprenta_rag.database.crud` to keep domain operations smaller and
more maintainable. The facade module re-exports these functions for backwards
compatibility.
"""

from __future__ import annotations

from typing import Any, List, Optional, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import Experiment
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_experiment(
    db: Session,
    name: str,
    experiment_type: Optional[str] = None,
    description: Optional[str] = None,
    disease: Optional[List[str]] = None,
    matrix: Optional[List[str]] = None,
    model_systems: Optional[List[str]] = None,
    targets: Optional[List[str]] = None,
    modality: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
    commit: bool = True,
) -> Experiment:
    """Create a new experiment in the database."""
    disease_list: List[str] = list(disease) if disease else []
    matrix_list: List[str] = list(matrix) if matrix else []
    model_systems_list: List[str] = list(model_systems) if model_systems else []
    targets_list: List[str] = list(targets) if targets else []
    modality_list: List[str] = list(modality) if modality else []
    experiment = Experiment(
        name=name,
        type=experiment_type,
        description=description,
        disease=cast(Any, disease_list),
        matrix=cast(Any, matrix_list),
        model_systems=cast(Any, model_systems_list),
        targets=cast(Any, targets_list),
        modality=cast(Any, modality_list),
        notion_page_id=notion_page_id,
    )

    db.add(experiment)

    if commit:
        db.commit()
        db.refresh(experiment)
        logger.info(
            "[CRUD][EXPERIMENT] Created experiment: %s (ID: %s)",
            experiment.name,
            experiment.id,
        )

    return experiment


def get_experiment_by_id(db: Session, experiment_id: UUID) -> Optional[Experiment]:
    """Get experiment by UUID."""
    return db.query(Experiment).filter(Experiment.id == experiment_id).first()


def get_experiment_by_name(db: Session, name: str) -> Optional[Experiment]:
    """Get experiment by name."""
    return db.query(Experiment).filter(Experiment.name == name).first()


def get_or_create_experiment(
    db: Session,
    name: str,
    **kwargs,
) -> tuple[Experiment, bool]:
    """
    Get existing experiment or create new one.

    Returns:
        Tuple of (Experiment, created: bool)
    """
    notion_page_id = kwargs.get("notion_page_id")

    # Try to find by notion_page_id first
    if notion_page_id:
        experiment = db.query(Experiment).filter(Experiment.notion_page_id == notion_page_id).first()
        if experiment:
            return experiment, False

    # Try to find by name
    experiment = get_experiment_by_name(db, name)
    if experiment:
        return experiment, False

    # Create new
    experiment = create_experiment(db, name=name, **kwargs)
    return experiment, True


