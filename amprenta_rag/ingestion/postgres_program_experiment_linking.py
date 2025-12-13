"""
Postgres program/experiment linking utilities.

Provides functions to link datasets to programs and experiments in Postgres,
including conversion between Notion page IDs and Postgres UUIDs.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import (
    Program,
    Experiment,
    program_dataset_assoc,
    experiment_dataset_assoc,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def find_program_by_notion_id(
    notion_page_id: str,
    db: Optional[Session] = None,
) -> Optional[Program]:
    """
    Find a Program in Postgres by Notion page ID.
    
    Args:
        notion_page_id: Notion page ID
        db: Optional database session
        
    Returns:
        Program model instance or None if not found
    """
    if db is None:
        with db_session() as db:
            return _find_program_by_notion_id_impl(notion_page_id, db)
    return _find_program_by_notion_id_impl(notion_page_id, db)


def find_experiment_by_notion_id(
    notion_page_id: str,
    db: Optional[Session] = None,
) -> Optional[Experiment]:
    """
    Find an Experiment in Postgres by Notion page ID.
    
    Args:
        notion_page_id: Notion page ID
        db: Optional database session
        
    Returns:
        Experiment model instance or None if not found
    """
    if db is None:
        with db_session() as db:
            return _find_experiment_by_notion_id_impl(notion_page_id, db)
    return _find_experiment_by_notion_id_impl(notion_page_id, db)


def convert_notion_ids_to_postgres_uuids(
    notion_program_ids: Optional[List[str]] = None,
    notion_experiment_ids: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> tuple[List[UUID], List[UUID]]:
    """
    Convert Notion page IDs to Postgres UUIDs for programs and experiments.
    
    Args:
        notion_program_ids: List of Notion program page IDs
        notion_experiment_ids: List of Notion experiment page IDs
        db: Optional database session
        
    Returns:
        Tuple of (program_uuids, experiment_uuids)
    """
    program_uuids: List[UUID] = []
    experiment_uuids: List[UUID] = []
    
    if db is None:
        with db_session() as db:
            return _convert_notion_ids_to_postgres_uuids_impl(
                notion_program_ids=notion_program_ids,
                notion_experiment_ids=notion_experiment_ids,
                db=db,
            )
    return _convert_notion_ids_to_postgres_uuids_impl(
        notion_program_ids=notion_program_ids,
        notion_experiment_ids=notion_experiment_ids,
        db=db,
    )


def link_dataset_to_programs_in_postgres(
    dataset_id: UUID,
    program_ids: List[UUID],
    db: Optional[Session] = None,
) -> int:
    """
    Link a dataset to programs in Postgres via association table.
    
    Args:
        dataset_id: Dataset UUID
        program_ids: List of Program UUIDs
        db: Optional database session
        
    Returns:
        Number of links created
    """
    if db is None:
        with db_session() as db:
            return _link_dataset_to_programs_in_postgres_impl(
                dataset_id=dataset_id,
                program_ids=program_ids,
                db=db,
            )
    return _link_dataset_to_programs_in_postgres_impl(
        dataset_id=dataset_id,
        program_ids=program_ids,
        db=db,
    )


def link_dataset_to_experiments_in_postgres(
    dataset_id: UUID,
    experiment_ids: List[UUID],
    db: Optional[Session] = None,
) -> int:
    """
    Link a dataset to experiments in Postgres via association table.
    
    Args:
        dataset_id: Dataset UUID
        experiment_ids: List of Experiment UUIDs
        db: Optional database session
        
    Returns:
        Number of links created
    """
    if db is None:
        with db_session() as db:
            return _link_dataset_to_experiments_in_postgres_impl(
                dataset_id=dataset_id,
                experiment_ids=experiment_ids,
                db=db,
            )
    return _link_dataset_to_experiments_in_postgres_impl(
        dataset_id=dataset_id,
        experiment_ids=experiment_ids,
        db=db,
    )


def _find_program_by_notion_id_impl(
    notion_page_id: str,
    db: Session,
) -> Optional[Program]:
    return (
        db.query(Program)
        .filter(Program.notion_page_id == notion_page_id)
        .first()
    )


def _find_experiment_by_notion_id_impl(
    notion_page_id: str,
    db: Session,
) -> Optional[Experiment]:
    return (
        db.query(Experiment)
        .filter(Experiment.notion_page_id == notion_page_id)
        .first()
    )


def _convert_notion_ids_to_postgres_uuids_impl(
    *,
    notion_program_ids: Optional[List[str]],
    notion_experiment_ids: Optional[List[str]],
    db: Session,
) -> tuple[List[UUID], List[UUID]]:
    program_uuids: List[UUID] = []
    experiment_uuids: List[UUID] = []

    if notion_program_ids:
        for notion_id in notion_program_ids:
            program = find_program_by_notion_id(notion_id, db=db)
            if program:
                program_uuids.append(program.id)
                logger.debug(
                    "[LINK][POSTGRES] Found Postgres Program %s for Notion ID %s",
                    program.id,
                    notion_id,
                )
            else:
                logger.warning(
                    "[LINK][POSTGRES] Program with Notion ID %s not found in Postgres",
                    notion_id,
                )

    if notion_experiment_ids:
        for notion_id in notion_experiment_ids:
            experiment = find_experiment_by_notion_id(notion_id, db=db)
            if experiment:
                experiment_uuids.append(experiment.id)
                logger.debug(
                    "[LINK][POSTGRES] Found Postgres Experiment %s for Notion ID %s",
                    experiment.id,
                    notion_id,
                )
            else:
                logger.warning(
                    "[LINK][POSTGRES] Experiment with Notion ID %s not found in Postgres",
                    notion_id,
                )

    return program_uuids, experiment_uuids


def _link_dataset_to_programs_in_postgres_impl(
    *,
    dataset_id: UUID,
    program_ids: List[UUID],
    db: Session,
) -> int:
    links_created = 0

    try:
        for program_id in program_ids:
            existing_link = db.execute(
                program_dataset_assoc.select().where(
                    program_dataset_assoc.c.program_id == program_id,
                    program_dataset_assoc.c.dataset_id == dataset_id,
                ),
            ).first()

            if existing_link:
                logger.debug(
                    "[LINK][POSTGRES] Dataset %s already linked to Program %s",
                    dataset_id,
                    program_id,
                )
                continue

            db.execute(
                program_dataset_assoc.insert().values(
                    program_id=program_id,
                    dataset_id=dataset_id,
                ),
            )
            links_created += 1

        if links_created > 0:
            db.commit()
            logger.info(
                "[LINK][POSTGRES] Linked dataset %s to %d programs",
                dataset_id,
                links_created,
            )

        return links_created

    except Exception as e:
        logger.error(
            "[LINK][POSTGRES] Error linking dataset %s to programs: %r",
            dataset_id,
            e,
        )
        db.rollback()
        raise


def _link_dataset_to_experiments_in_postgres_impl(
    *,
    dataset_id: UUID,
    experiment_ids: List[UUID],
    db: Session,
) -> int:
    links_created = 0

    try:
        for experiment_id in experiment_ids:
            existing_link = db.execute(
                experiment_dataset_assoc.select().where(
                    experiment_dataset_assoc.c.experiment_id == experiment_id,
                    experiment_dataset_assoc.c.dataset_id == dataset_id,
                ),
            ).first()

            if existing_link:
                logger.debug(
                    "[LINK][POSTGRES] Dataset %s already linked to Experiment %s",
                    dataset_id,
                    experiment_id,
                )
                continue

            db.execute(
                experiment_dataset_assoc.insert().values(
                    experiment_id=experiment_id,
                    dataset_id=dataset_id,
                ),
            )
            links_created += 1

        if links_created > 0:
            db.commit()
            logger.info(
                "[LINK][POSTGRES] Linked dataset %s to %d experiments",
                dataset_id,
                links_created,
            )

        return links_created

    except Exception as e:
        logger.error(
            "[LINK][POSTGRES] Error linking dataset %s to experiments: %r",
            dataset_id,
            e,
        )
        db.rollback()
        raise


def link_dataset_to_programs_and_experiments_in_postgres(
    dataset_id: UUID,
    program_ids: Optional[List[UUID]] = None,
    experiment_ids: Optional[List[UUID]] = None,
    notion_program_ids: Optional[List[str]] = None,
    notion_experiment_ids: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> dict[str, int]:
    """
    Link a dataset to programs and experiments in Postgres.
    
    Accepts either Postgres UUIDs directly or Notion page IDs (which will be converted).
    
    Args:
        dataset_id: Dataset UUID
        program_ids: List of Program UUIDs (direct)
        experiment_ids: List of Experiment UUIDs (direct)
        notion_program_ids: List of Notion program page IDs (will be converted)
        notion_experiment_ids: List of Notion experiment page IDs (will be converted)
        db: Optional database session
        
    Returns:
        Dict with counts: {'programs_linked': int, 'experiments_linked': int}
    """
    # Convert Notion IDs if provided
    if notion_program_ids or notion_experiment_ids:
        pg_program_ids, pg_experiment_ids = convert_notion_ids_to_postgres_uuids(
            notion_program_ids=notion_program_ids,
            notion_experiment_ids=notion_experiment_ids,
            db=db,
        )
        
        # Merge with direct UUIDs
        if program_ids:
            pg_program_ids = list(set(pg_program_ids + program_ids))
        if experiment_ids:
            pg_experiment_ids = list(set(pg_experiment_ids + experiment_ids))
    else:
        pg_program_ids = program_ids or []
        pg_experiment_ids = experiment_ids or []
    
    results = {
        "programs_linked": 0,
        "experiments_linked": 0,
    }
    
    # Link programs
    if pg_program_ids:
        results["programs_linked"] = link_dataset_to_programs_in_postgres(
            dataset_id=dataset_id,
            program_ids=pg_program_ids,
            db=db,
        )
    
    # Link experiments
    if pg_experiment_ids:
        results["experiments_linked"] = link_dataset_to_experiments_in_postgres(
            dataset_id=dataset_id,
            experiment_ids=pg_experiment_ids,
            db=db,
        )
    
    return results

