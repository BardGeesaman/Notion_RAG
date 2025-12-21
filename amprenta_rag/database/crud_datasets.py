"""
Dataset (and related Feature) CRUD operations.

Split out of `amprenta_rag.database.crud` to keep domain operations smaller and
more maintainable. The facade module re-exports these functions for backwards
compatibility.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from sqlalchemy import and_, or_
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Feature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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
    **kwargs,
) -> tuple[Dataset, bool]:
    """
    Get existing dataset or create new one.

    Returns:
        Tuple of (Dataset, created: bool)
    """
    notion_page_id = kwargs.get("notion_page_id")

    # Try to find by notion_page_id first
    if notion_page_id:
        dataset = db.query(Dataset).filter(Dataset.notion_page_id == notion_page_id).first()
        if dataset:
            return dataset, False

    # Try to find by name and omics_type
    dataset = db.query(Dataset).filter(and_(Dataset.name == name, Dataset.omics_type == omics_type)).first()
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
# FEATURE CRUD (owned by the datasets domain in this split)
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
    feature_type: str,
) -> Optional[Feature]:
    """Get feature by name and type."""
    return (
        db.query(Feature)
        .filter(
            and_(
                or_(
                    Feature.name == name,
                    Feature.normalized_name == name,
                ),
                Feature.feature_type == feature_type,
            )
        )
        .first()
    )


def get_or_create_feature(
    db: Session,
    name: str,
    feature_type: str,
    normalized_name: Optional[str] = None,
    **kwargs,
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
        **kwargs,
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
    features: List[Feature] = []

    for i, data in enumerate(features_data):
        feature, _created = get_or_create_feature(
            db,
            name=data["name"],
            feature_type=data["feature_type"],
            normalized_name=data.get("normalized_name"),
            commit=False,
        )
        features.append(feature)

        # Commit in batches
        if (i + 1) % batch_size == 0:
            db.commit()
            logger.debug(
                "[CRUD][FEATURE] Committed batch of %d features",
                batch_size,
            )

    # Final commit
    db.commit()
    logger.info(
        "[CRUD][FEATURE] Created/retrieved %d features",
        len(features),
    )

    return features


