"""
Postgres integration for ingestion pipelines.

This module provides functions for ingestion modules to write directly to
Postgres as the system of record, replacing Notion as the primary database.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.session import db_session
from amprenta_rag.database.crud import (
    get_or_create_dataset,
    link_dataset_to_features,
    bulk_get_or_create_features,
)
from amprenta_rag.database.models import Dataset, Feature
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_or_update_dataset_in_postgres(
    name: str,
    omics_type: OmicsType,
    file_paths: Optional[List[str]] = None,
    description: Optional[str] = None,
    organism: Optional[List[str]] = None,
    sample_type: Optional[List[str]] = None,
    disease: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
    external_ids: Optional[Dict[str, str]] = None,
    auto_ingest: bool = False,
    db: Optional[Session] = None,
) -> Dataset:
    """
    Create or update a dataset in Postgres.

    Args:
        name: Dataset name
        omics_type: OmicsType enum value
        file_paths: List of file paths
        description: Dataset description
        organism: List of organism names
        sample_type: List of sample types
        disease: List of disease terms
        notion_page_id: Notion page ID for migration compatibility
        program_ids: List of program IDs (Notion or Postgres UUIDs)
        experiment_ids: List of experiment IDs (Notion or Postgres UUIDs)
        external_ids: Dict of external IDs (e.g., {"mw_study_id": "ST001234"})
        auto_ingest: Whether to trigger automatic ingestion (not yet implemented)
        db: Database session (optional, will create if not provided)

    Returns:
        Dataset object
    """
    if db is None:
        with db_session() as db:
            return _create_or_update_dataset_in_postgres_impl(
                name=name,
                omics_type=omics_type,
                file_paths=file_paths,
                description=description,
                organism=organism,
                sample_type=sample_type,
                disease=disease,
                notion_page_id=notion_page_id,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
                external_ids=external_ids,
                auto_ingest=auto_ingest,
                db=db,
            )
    return _create_or_update_dataset_in_postgres_impl(
        name=name,
        omics_type=omics_type,
        file_paths=file_paths,
        description=description,
        organism=organism,
        sample_type=sample_type,
        disease=disease,
        notion_page_id=notion_page_id,
        program_ids=program_ids,
        experiment_ids=experiment_ids,
        external_ids=external_ids,
        auto_ingest=auto_ingest,
        db=db,
    )


def link_features_to_dataset_postgres(
    dataset_id: UUID,
    feature_names: List[str],
    feature_type: str,
    normalized_names: Optional[List[str]] = None,
    db: Optional[Session] = None,
) -> List[Feature]:
    """
    Link features to a dataset in Postgres.

    Creates features if they don't exist and links them to the dataset.

    Args:
        dataset_id: Dataset UUID
        feature_names: List of feature names
        feature_type: Feature type (gene, protein, metabolite, lipid)
        normalized_names: Optional list of normalized names (parallel to feature_names)
        db: Database session (optional)

    Returns:
        List of Feature objects
    """
    if db is None:
        with db_session() as db:
            return _link_features_to_dataset_postgres_impl(
                dataset_id=dataset_id,
                feature_names=feature_names,
                feature_type=feature_type,
                normalized_names=normalized_names,
                db=db,
            )
    return _link_features_to_dataset_postgres_impl(
        dataset_id=dataset_id,
        feature_names=feature_names,
        feature_type=feature_type,
        normalized_names=normalized_names,
        db=db,
    )


def get_postgres_session() -> Session:
    """Get a Postgres database session."""
    return db_session().__enter__()


def embed_dataset_with_postgres_metadata(
    dataset_id: UUID,
    dataset_name: str,
    species_or_features: List[str],
    omics_type: OmicsType,
    signature_matches: Optional[List[Any]] = None,
    notion_page_id: Optional[str] = None,
) -> None:
    """
    Embed a dataset into Pinecone using Postgres metadata.

    This is the Postgres-aware counterpart to the Notion-only embedding functions
    (e.g., `embed_lipidomics_dataset`, `embed_metabolomics_dataset`). It is
    intentionally implemented as a thin shim so ingestion pipelines can depend
    on a stable API while we finalize the Postgres-centric embedding design.

    Current behavior:
    - Logs that the function is not yet fully implemented.
    - Raises NotImplementedError so ingestion callers can fall back to their
      existing Notion-based embedding inside their try/except blocks.

    Args:
        dataset_id: Postgres dataset UUID.
        dataset_name: Human-readable dataset name.
        species_or_features: List of feature identifiers (e.g., metabolites, lipids).
        omics_type: OmicsType enum describing the dataset.
        signature_matches: Optional list of signature match results.
        notion_page_id: Optional Notion page ID linked to this dataset.
    """
    logger.info(
        "[POSTGRES][EMBED] embed_dataset_with_postgres_metadata() called for dataset %s (ID=%s, omics_type=%s); "
        "Postgres-aware embedding is not yet implemented, callers should fall back to Notion-based embedding.",
        dataset_name,
        dataset_id,
        getattr(omics_type, "value", str(omics_type)),
    )
    raise NotImplementedError(
        "embed_dataset_with_postgres_metadata is not yet implemented; "
        "ingestion pipelines should catch this and fall back to Notion-based embedding."
    )


def _create_or_update_dataset_in_postgres_impl(
    *,
    name: str,
    omics_type: OmicsType,
    file_paths: Optional[List[str]],
    description: Optional[str],
    organism: Optional[List[str]],
    sample_type: Optional[List[str]],
    disease: Optional[List[str]],
    notion_page_id: Optional[str],
    program_ids: Optional[List[str]],
    experiment_ids: Optional[List[str]],
    external_ids: Optional[Dict[str, str]],
    auto_ingest: bool,
    db: Session,
) -> Dataset:
    omics_type_str = omics_type.value if hasattr(omics_type, "value") else str(omics_type)

    dataset, created = get_or_create_dataset(
        db,
        name=name,
        omics_type=omics_type_str,
        description=description,
        file_paths=file_paths or [],
        organism=organism or [],
        sample_type=sample_type or [],
        disease=disease or [],
        notion_page_id=notion_page_id,
    )

    if external_ids:
        if dataset.external_ids:
            dataset.external_ids.update(external_ids)
        else:
            dataset.external_ids = external_ids
        db.commit()

    action = "Created" if created else "Found existing"
    logger.info(
        "[POSTGRES][DATASET] %s dataset: %s (ID: %s, external_ids: %s)",
        action,
        dataset.name,
        dataset.id,
        external_ids,
    )

    return dataset


def _link_features_to_dataset_postgres_impl(
    *,
    dataset_id: UUID,
    feature_names: List[str],
    feature_type: str,
    normalized_names: Optional[List[str]],
    db: Session,
) -> List[Feature]:
    features_data = []
    for i, name in enumerate(feature_names):
        norm_name = normalized_names[i] if normalized_names and i < len(normalized_names) else name
        features_data.append(
            {
                "name": name,
                "feature_type": feature_type,
                "normalized_name": norm_name,
            }
        )

    features = bulk_get_or_create_features(db, features_data)

    feature_ids: List[UUID] = [cast(UUID, f.id) for f in features if f.id is not None]
    link_dataset_to_features(db, dataset_id, feature_ids)

    logger.info(
        "[POSTGRES][FEATURES] Linked %d %s features to dataset %s",
        len(features),
        feature_type,
        dataset_id,
    )

    return features
