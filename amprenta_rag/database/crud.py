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

# ==============================================================================
# PROGRAM CRUD
# ==============================================================================
# ==============================================================================
# EXPERIMENT CRUD
# ==============================================================================
# ==============================================================================
# DATASET CRUD
# ==============================================================================

def create_dataset(
    db: Session,
    name: str,
    omics_type: str,
    description: Optional[str] = None,
    file_paths: Optional[List[str]] = None,
    file_urls: Optional[List[str]] = None,
    organism: Optional[List[str]] = None,
    sample_type: Optional[List[str]] = None,
    disease: Optional[List[str]] = None,
    methods: Optional[str] = None,
    summary: Optional[str] = None,
    results: Optional[str] = None,
    conclusions: Optional[str] = None,
    dataset_source_type: Optional[str] = None,
    data_origin: Optional[str] = None,
    notion_page_id: Optional[str] = None,
    external_ids: Optional[Dict[str, Any]] = None,
    commit: bool = True,
) -> Dataset:
    """Create a new dataset in the database."""
    file_paths_list: List[str] = list(file_paths) if file_paths else []
    file_urls_list: List[str] = list(file_urls) if file_urls else []
    organism_list: List[str] = list(organism) if organism else []
    sample_type_list: List[str] = list(sample_type) if sample_type else []
    disease_list: List[str] = list(disease) if disease else []
    dataset = Dataset(
        name=name,
        omics_type=omics_type,
        description=description,
        file_paths=cast(Any, file_paths_list),
        file_urls=cast(Any, file_urls_list),
        organism=cast(Any, organism_list),
        sample_type=cast(Any, sample_type_list),
        disease=cast(Any, disease_list),
        methods=methods,
        summary=summary,
        results=results,
        conclusions=conclusions,
        dataset_source_type=dataset_source_type,
        data_origin=data_origin,
        notion_page_id=notion_page_id,
        external_ids=external_ids or {},
    )

    db.add(dataset)

    if commit:
        db.commit()
        db.refresh(dataset)
        logger.info(
            "[CRUD][DATASET] Created dataset: %s (ID: %s, Type: %s)",
            dataset.name,
            dataset.id,
            dataset.omics_type,
        )

    return dataset


def get_dataset_by_id(db: Session, dataset_id: UUID) -> Optional[Dataset]:
    """Get dataset by UUID."""
    return db.query(Dataset).filter(Dataset.id == dataset_id).first()


def get_dataset_by_name(db: Session, name: str) -> Optional[Dataset]:
    """Get dataset by name."""
    return db.query(Dataset).filter(Dataset.name == name).first()


def get_or_create_dataset(
    db: Session,
    name: str,
    omics_type: str,
    **kwargs
) -> tuple[Dataset, bool]:
    """
    Get existing dataset or create new one.

    Returns:
        Tuple of (Dataset, created: bool)
    """
    notion_page_id = kwargs.get('notion_page_id')

    # Try to find by notion_page_id first
    if notion_page_id:
        dataset = db.query(Dataset).filter(
            Dataset.notion_page_id == notion_page_id
        ).first()
        if dataset:
            return dataset, False

    # Try to find by name and omics_type
    dataset = db.query(Dataset).filter(
        and_(Dataset.name == name, Dataset.omics_type == omics_type)
    ).first()
    if dataset:
        return dataset, False

    # Create new
    dataset = create_dataset(db, name=name, omics_type=omics_type, **kwargs)
    return dataset, True


def list_datasets(
    db: Session,
    omics_type: Optional[str] = None,
    limit: Optional[int] = None,
    offset: int = 0,
) -> List[Dataset]:
    """List datasets with optional filtering."""
    query = db.query(Dataset)

    if omics_type:
        query = query.filter(Dataset.omics_type == omics_type)

    query = query.order_by(Dataset.name)

    if limit:
        query = query.limit(limit).offset(offset)

    return query.all()


# ==============================================================================
# FEATURE CRUD
# ==============================================================================

def create_feature(
    db: Session,
    name: str,
    feature_type: str,
    normalized_name: Optional[str] = None,
    aliases: Optional[List[str]] = None,
    notion_page_id: Optional[str] = None,
    external_ids: Optional[Dict[str, Any]] = None,
    commit: bool = True,
) -> Feature:
    """Create a new feature in the database."""
    aliases_list: List[str] = list(aliases) if aliases else []
    feature = Feature(
        name=name,
        feature_type=feature_type,
        normalized_name=normalized_name or name,
        aliases=cast(Any, aliases_list),
        notion_page_id=notion_page_id,
        external_ids=external_ids or {},
    )

    db.add(feature)

    if commit:
        db.commit()
        db.refresh(feature)
        logger.debug(
            "[CRUD][FEATURE] Created feature: %s (Type: %s, ID: %s)",
            feature.name,
            feature.feature_type,
            feature.id,
        )

    return feature


def get_feature_by_id(db: Session, feature_id: UUID) -> Optional[Feature]:
    """Get feature by UUID."""
    return db.query(Feature).filter(Feature.id == feature_id).first()


def get_feature_by_name_and_type(
    db: Session,
    name: str,
    feature_type: str
) -> Optional[Feature]:
    """Get feature by name and type."""
    return db.query(Feature).filter(
        and_(
            or_(
                Feature.name == name,
                Feature.normalized_name == name
            ),
            Feature.feature_type == feature_type
        )
    ).first()


def get_or_create_feature(
    db: Session,
    name: str,
    feature_type: str,
    normalized_name: Optional[str] = None,
    **kwargs
) -> tuple[Feature, bool]:
    """
    Get existing feature or create new one.

    Uses normalized_name for matching if provided.

    Returns:
        Tuple of (Feature, created: bool)
    """
    # Search by normalized name first
    search_name = normalized_name or name

    feature = get_feature_by_name_and_type(db, search_name, feature_type)
    if feature:
        return feature, False

    # Create new
    feature = create_feature(
        db,
        name=name,
        feature_type=feature_type,
        normalized_name=normalized_name,
        **kwargs
    )
    return feature, True


def bulk_get_or_create_features(
    db: Session,
    features_data: List[Dict[str, Any]],
    batch_size: int = 100,
) -> List[Feature]:
    """
    Bulk get or create features for better performance.

    Args:
        db: Database session
        features_data: List of dicts with keys: name, feature_type, normalized_name (optional)
        batch_size: Batch size for committing

    Returns:
        List of Feature objects
    """
    features = []

    for i, data in enumerate(features_data):
        feature, created = get_or_create_feature(
            db,
            name=data['name'],
            feature_type=data['feature_type'],
            normalized_name=data.get('normalized_name'),
            commit=False,
        )
        features.append(feature)

        # Commit in batches
        if (i + 1) % batch_size == 0:
            db.commit()
            logger.debug(
                "[CRUD][FEATURE] Committed batch of %d features",
                batch_size
            )

    # Final commit
    db.commit()
    logger.info(
        "[CRUD][FEATURE] Created/retrieved %d features",
        len(features)
    )

    return features


# ==============================================================================
# SIGNATURE CRUD
# ==============================================================================

def create_signature(
    db: Session,
    name: str,
    description: Optional[str] = None,
    modalities: Optional[List[str]] = None,
    short_id: Optional[str] = None,
    biomarker_role: Optional[List[str]] = None,
    phenotype_axes: Optional[List[str]] = None,
    data_ownership: Optional[str] = None,
    notion_page_id: Optional[str] = None,
    commit: bool = True,
) -> Signature:
    """Create a new signature in the database."""
    modalities_list: List[str] = list(modalities) if modalities else []
    biomarker_list: List[str] = list(biomarker_role) if biomarker_role else []
    phenotype_list: List[str] = list(phenotype_axes) if phenotype_axes else []
    signature = Signature(
        name=name,
        description=description,
        modalities=cast(Any, modalities_list),
        short_id=short_id,
        biomarker_role=cast(Any, biomarker_list),
        phenotype_axes=cast(Any, phenotype_list),
        data_ownership=data_ownership,
        notion_page_id=notion_page_id,
    )

    db.add(signature)

    if commit:
        db.commit()
        db.refresh(signature)
        logger.info(
            "[CRUD][SIGNATURE] Created signature: %s (ID: %s)",
            signature.name,
            signature.id,
        )

    return signature


def create_signature_component(
    db: Session,
    signature_id: UUID,
    feature_id: Optional[UUID],
    feature_name: str,
    feature_type: str,
    direction: Optional[str] = None,
    weight: float = 1.0,
    commit: bool = True,
) -> SignatureComponent:
    """Create a signature component."""
    weight_value: float | None = float(weight) if weight is not None else None
    component = SignatureComponent(
        signature_id=cast(Any, signature_id),
        feature_id=cast(Any, feature_id),
        feature_name=feature_name,
        feature_type=feature_type,
        direction=direction,
        weight=cast(Any, weight_value),
    )

    db.add(component)

    if commit:
        db.commit()
        db.refresh(component)

    return component


def get_signature_by_id(db: Session, signature_id: UUID) -> Optional[Signature]:
    """Get signature by UUID."""
    return db.query(Signature).filter(Signature.id == signature_id).first()


def get_signature_by_name(db: Session, name: str) -> Optional[Signature]:
    """Get signature by name."""
    return db.query(Signature).filter(Signature.name == name).first()


def get_or_create_signature(
    db: Session,
    name: str,
    **kwargs
) -> tuple[Signature, bool]:
    """
    Get existing signature or create new one.

    Returns:
        Tuple of (Signature, created: bool)
    """
    notion_page_id = kwargs.get('notion_page_id')

    # Try to find by notion_page_id first
    if notion_page_id:
        signature = db.query(Signature).filter(
            Signature.notion_page_id == notion_page_id
        ).first()
        if signature:
            return signature, False

    # Try to find by name
    signature = get_signature_by_name(db, name)
    if signature:
        return signature, False

    # Create new
    signature = create_signature(db, name=name, **kwargs)
    return signature, True


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


def link_dataset_to_signature(
    db: Session,
    dataset_id: UUID,
    signature_id: UUID,
    match_score: Optional[float] = None,
    commit: bool = True,
) -> None:
    """Link a dataset to a signature with an optional match score."""
    dataset = get_dataset_by_id(db, dataset_id)
    signature = get_signature_by_id(db, signature_id)

    if not dataset or not signature:
        logger.error(
            "[CRUD][LINK] Cannot link dataset %s to signature %s: not found",
            dataset_id,
            signature_id,
        )
        return

    if signature not in dataset.signatures:
        dataset.signatures.append(signature)

        # If match_score is provided, update the junction table
        # This requires direct SQL manipulation - simplified version for now
        if match_score is not None and hasattr(dataset, "signature_match_score"):
            setattr(dataset, "signature_match_score", float(match_score))

        if commit:
            db.commit()
            logger.debug(
                "[CRUD][LINK] Linked dataset %s to signature %s (score: %.3f)",
                dataset.name,
                signature.name,
                match_score or 0.0,
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

