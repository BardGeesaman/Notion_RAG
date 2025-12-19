"""
Postgres-based signature linking.

Handles linking signatures to datasets in Postgres using the dataset_signature_assoc table.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import (
    Dataset as DatasetModel,
    Signature as SignatureModel,
    dataset_signature_assoc,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def link_signature_to_dataset_in_postgres(
    signature_id: UUID,
    dataset_id: UUID,
    match_score: Optional[float] = None,
    db: Optional[Session] = None,
) -> bool:
    """
    Link a signature to a dataset in Postgres.

    Creates an entry in the dataset_signature_assoc table with optional match_score.

    Args:
        signature_id: Postgres UUID of the signature
        dataset_id: Postgres UUID of the dataset
        match_score: Optional match score (overlap fraction, etc.)
        db: Optional database session (will create one if not provided)

    Returns:
        True if linking succeeded, False otherwise
    """
    if db is None:
        with db_session() as db:
            return _link_signature_to_dataset_impl(
                signature_id=signature_id,
                dataset_id=dataset_id,
                match_score=match_score,
                db=db,
            )
    return _link_signature_to_dataset_impl(
        signature_id=signature_id,
        dataset_id=dataset_id,
        match_score=match_score,
        db=db,
    )


def get_dataset_signatures_from_postgres(
    dataset_id: UUID,
) -> List[SignatureModel]:
    """
    Get all signatures linked to a dataset from Postgres.

    Args:
        dataset_id: Postgres UUID of the dataset

    Returns:
        List of SignatureModel instances linked to the dataset
    """
    with db_session() as db:
        try:
            dataset = (
                db.query(DatasetModel)
                .filter(DatasetModel.id == dataset_id)
                .first()
            )

            if not dataset:
                logger.warning(
                    "[SIGNATURE-LINK] Dataset %s not found in Postgres",
                    dataset_id,
                )
                return []

            signatures = dataset.signatures or []

            logger.debug(
                "[SIGNATURE-LINK] Found %d signature(s) linked to dataset %s",
                len(signatures),
                dataset_id,
            )

            return signatures

        except Exception as e:
            logger.error(
                "[SIGNATURE-LINK] Error getting signatures for dataset %s: %r",
                dataset_id,
                e,
            )
            return []


def get_signature_datasets_from_postgres(
    signature_id: UUID,
) -> List[DatasetModel]:
    """
    Get all datasets linked to a signature from Postgres.

    Args:
        signature_id: Postgres UUID of the signature

    Returns:
        List of DatasetModel instances linked to the signature
    """
    with db_session() as db:
        try:
            signature = (
                db.query(SignatureModel)
                .filter(SignatureModel.id == signature_id)
                .first()
            )

            if not signature:
                logger.warning(
                    "[SIGNATURE-LINK] Signature %s not found in Postgres",
                    signature_id,
                )
                return []

            datasets = signature.datasets or []

            logger.debug(
                "[SIGNATURE-LINK] Found %d dataset(s) linked to signature %s",
                len(datasets),
                signature_id,
            )

            return datasets

        except Exception as e:
            logger.error(
                "[SIGNATURE-LINK] Error getting datasets for signature %s: %r",
                signature_id,
                e,
            )
            return []


def unlink_signature_from_dataset(
    signature_id: UUID,
    dataset_id: UUID,
    db: Optional[Session] = None,
) -> bool:
    """
    Unlink a signature from a dataset in Postgres.

    Args:
        signature_id: Postgres UUID of the signature
        dataset_id: Postgres UUID of the dataset
        db: Optional database session

    Returns:
        True if unlinking succeeded, False otherwise
    """
    if db is None:
        with db_session() as db:
            return _unlink_signature_from_dataset_impl(
                signature_id=signature_id,
                dataset_id=dataset_id,
                db=db,
            )
    return _unlink_signature_from_dataset_impl(
        signature_id=signature_id,
        dataset_id=dataset_id,
        db=db,
    )


def _link_signature_to_dataset_impl(
    *,
    signature_id: UUID,
    dataset_id: UUID,
    match_score: Optional[float],
    db: Session,
) -> bool:
    try:
        existing = (
            db.query(dataset_signature_assoc)
            .filter(
                dataset_signature_assoc.c.dataset_id == dataset_id,
                dataset_signature_assoc.c.signature_id == signature_id,
            )
            .first()
        )

        if existing:
            if match_score is not None:
                db.execute(
                    dataset_signature_assoc.update()
                    .where(
                        dataset_signature_assoc.c.dataset_id == dataset_id,
                        dataset_signature_assoc.c.signature_id == signature_id,
                    )
                    .values(match_score=match_score)
                )
                db.commit()
                logger.debug(
                    "[SIGNATURE-LINK] Updated match score for signature %s -> dataset %s: %.3f",
                    signature_id,
                    dataset_id,
                    match_score,
                )
            return True

        insert_stmt = dataset_signature_assoc.insert().values(
            dataset_id=dataset_id,
            signature_id=signature_id,
            match_score=match_score,
        )
        db.execute(insert_stmt)
        db.commit()

        logger.info(
            "[SIGNATURE-LINK] Linked signature %s to dataset %s (score: %s)",
            signature_id,
            dataset_id,
            match_score if match_score is not None else "N/A",
        )

        return True

    except Exception as e:
        logger.error(
            "[SIGNATURE-LINK] Error linking signature %s to dataset %s: %r",
            signature_id,
            dataset_id,
            e,
        )
        db.rollback()
        return False


def _unlink_signature_from_dataset_impl(
    *,
    signature_id: UUID,
    dataset_id: UUID,
    db: Session,
) -> bool:
    try:
        delete_stmt = dataset_signature_assoc.delete().where(
            dataset_signature_assoc.c.dataset_id == dataset_id,
            dataset_signature_assoc.c.signature_id == signature_id,
        )
        result = db.execute(delete_stmt)
        db.commit()

        if result.rowcount > 0:
            logger.info(
                "[SIGNATURE-LINK] Unlinked signature %s from dataset %s",
                signature_id,
                dataset_id,
            )
            return True

        logger.debug(
            "[SIGNATURE-LINK] No association found between signature %s and dataset %s",
            signature_id,
            dataset_id,
        )
        return False

    except Exception as e:
        logger.error(
            "[SIGNATURE-LINK] Error unlinking signature %s from dataset %s: %r",
            signature_id,
            dataset_id,
            e,
        )
        db.rollback()
        return False

