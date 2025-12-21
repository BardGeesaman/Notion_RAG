"""
CRUD operations for Postgres database.

This module provides a repository pattern implementation for all database
operations. All ingestion modules should use these functions instead of
direct Notion API calls or raw SQLAlchemy queries.

Key Principles:
- Idempotent operations (safe to run multiple times)
- Clear error handling
- Comprehensive logging
- UUID-based IDs
- Maintains notion_page_id for migration compatibility
"""

from __future__ import annotations

from typing import List, Optional, Dict, Any, Sequence, cast
from uuid import UUID
import uuid

from sqlalchemy import and_, or_
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    Program,
    Experiment,
    Dataset,
    Feature,
    Signature,
    SignatureComponent,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Domain splits (facade re-exports keep backwards compatibility)
from .crud_programs import (  # noqa: E402
    create_program,
    get_program_by_id,
    get_program_by_name,
    get_program_by_notion_id,
    get_or_create_program,
    list_programs,
)
from .crud_experiments import (  # noqa: E402
    create_experiment,
    get_experiment_by_id,
    get_experiment_by_name,
    get_or_create_experiment,
)
from .crud_datasets import (  # noqa: E402
    create_dataset,
    get_dataset_by_id,
    get_dataset_by_name,
    get_or_create_dataset,
    list_datasets,
    create_feature,
    get_feature_by_id,
    get_feature_by_name_and_type,
    get_or_create_feature,
    bulk_get_or_create_features,
)
from .crud_signatures import (  # noqa: E402
    create_signature,
    create_signature_component,
    get_signature_by_id,
    get_signature_by_name,
    get_or_create_signature,
    link_dataset_to_signature,
)

# ==============================================================================
# PROGRAM CRUD
# ==============================================================================
# ==============================================================================
# EXPERIMENT CRUD
# ==============================================================================
# ==============================================================================
# DATASET CRUD
# ==============================================================================
# ==============================================================================
# SIGNATURE CRUD
# ==============================================================================
# ==============================================================================
# RELATIONSHIP MANAGEMENT
# ==============================================================================

def link_dataset_to_program(
    db: Session,
    dataset_id: UUID,
    program_id: UUID,
    commit: bool = True,
) -> None:
    """Link a dataset to a program."""
    dataset = get_dataset_by_id(db, dataset_id)
    program = get_program_by_id(db, program_id)

    if not dataset or not program:
        logger.error(
            "[CRUD][LINK] Cannot link dataset %s to program %s: not found",
            dataset_id,
            program_id,
        )
        return

    if program not in dataset.programs:
        dataset.programs.append(program)

        if commit:
            db.commit()
            logger.debug(
                "[CRUD][LINK] Linked dataset %s to program %s",
                dataset.name,
                program.name,
            )


def link_dataset_to_experiment(
    db: Session,
    dataset_id: UUID,
    experiment_id: UUID,
    commit: bool = True,
) -> None:
    """Link a dataset to an experiment."""
    dataset = get_dataset_by_id(db, dataset_id)
    experiment = get_experiment_by_id(db, experiment_id)

    if not dataset or not experiment:
        logger.error(
            "[CRUD][LINK] Cannot link dataset %s to experiment %s: not found",
            dataset_id,
            experiment_id,
        )
        return

    if experiment not in dataset.experiments:
        dataset.experiments.append(experiment)

        if commit:
            db.commit()
            logger.debug(
                "[CRUD][LINK] Linked dataset %s to experiment %s",
                dataset.name,
                experiment.name,
            )


def link_dataset_to_features(
    db: Session,
    dataset_id: UUID,
    feature_ids: List[UUID],
    commit: bool = True,
) -> None:
    """Link a dataset to multiple features."""
    dataset = get_dataset_by_id(db, dataset_id)

    if not dataset:
        logger.error("[CRUD][LINK] Dataset %s not found", dataset_id)
        return

    for feature_id in feature_ids:
        feature = get_feature_by_id(db, feature_id)
        if feature and feature not in dataset.features:
            dataset.features.append(feature)

    if commit:
        db.commit()
        logger.info(
            "[CRUD][LINK] Linked %d features to dataset %s",
            len(feature_ids),
            dataset.name,
        )


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def get_database_stats(db: Session) -> Dict[str, int]:
    """Get counts of all major entities in the database."""
    return {
        'programs': db.query(Program).count(),
        'experiments': db.query(Experiment).count(),
        'datasets': db.query(Dataset).count(),
        'features': db.query(Feature).count(),
        'signatures': db.query(Signature).count(),
        'signature_components': db.query(SignatureComponent).count(),
    }

