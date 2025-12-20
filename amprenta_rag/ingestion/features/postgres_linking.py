"""
Postgres-based feature linking.

Links metabolite features to Postgres entities (Literature, Email, Dataset, Experiment).
This is the Postgres-first alternative to link_features_to_notion_items().
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.database.models import Email
from amprenta_rag.database.models import Feature as FeatureModel
from amprenta_rag.database.models import Literature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def link_features_to_postgres_items(
    feature_names: List[str],
    item_id: UUID,
    item_type: str,
    db: Optional[Session] = None,
) -> None:
    """
    Link metabolite features to a Postgres entity (dataset, literature, email, experiment).

    For each feature:
    1. Finds or creates the corresponding Feature record in Postgres
    2. Links the feature to the target entity:
       - For datasets: Uses dataset_feature_assoc relationship
       - For literature/email: Stores feature mentions in semantic_metadata JSON

    Args:
        feature_names: List of normalized metabolite names
        item_id: UUID of the target entity in Postgres
        item_type: One of "dataset", "literature", "email", "experiment"
    """
    if db is None:
        db_gen = get_db()
        session = next(db_gen)
        try:
            return link_features_to_postgres_items(feature_names, item_id, item_type, db=session)
        finally:
            session.close()
            next(db_gen, None)

    if not feature_names:
        logger.debug(
            "[INGEST][FEATURES] No features to link for %s %s",
            item_type,
            item_id,
        )
        return

    logger.info(
        "[INGEST][FEATURES] Linking %d feature(s) to %s %s",
        len(feature_names),
        item_type,
        item_id,
    )

    linked_count = 0
    feature_mentions: List[str] = []

    for feature_name in feature_names:
        try:
            # Find or create feature
            feature = _find_or_create_feature_impl(feature_name, "metabolite", db)
            if not feature:
                logger.warning(
                    "[INGEST][FEATURES] Could not create/find feature %s",
                    feature_name,
                )
                continue

            # Link based on item_type
            if item_type == "dataset":
                # Use dataset_feature_assoc relationship
                dataset = db.query(DatasetModel).filter(DatasetModel.id == item_id).first()
                if dataset and feature not in dataset.features:
                    dataset.features.append(feature)
                    linked_count += 1
                    logger.debug(
                        "[INGEST][FEATURES] Linked feature %s to dataset %s",
                        feature_name,
                        item_id,
                    )
            elif item_type in ("literature", "email"):
                # Store feature mentions in semantic_metadata
                feature_mentions.append(feature_name)
                linked_count += 1
                logger.debug(
                    "[INGEST][FEATURES] Added feature mention %s to %s %s",
                    feature_name,
                    item_type,
                    item_id,
                )
            elif item_type == "experiment":
                # For experiments, we can link through datasets
                # For now, store in metadata similar to literature/email
                feature_mentions.append(feature_name)
                linked_count += 1
                logger.debug(
                    "[INGEST][FEATURES] Added feature mention %s to experiment %s",
                    feature_name,
                    item_id,
                )
            else:
                logger.warning(
                    "[INGEST][FEATURES] Unknown item_type %s for %s",
                    item_type,
                    item_id,
                )

        except Exception as e:
            logger.warning(
                "[INGEST][FEATURES] Error linking feature %s to %s %s: %r",
                feature_name,
                item_type,
                item_id,
                e,
            )

    # Update semantic_metadata for literature/email/experiment
    if feature_mentions and item_type in ("literature", "email", "experiment"):
        _update_semantic_metadata_with_features(item_id, item_type, feature_mentions, db=db)

    logger.info(
        "[INGEST][FEATURES] Successfully linked %d/%d features to %s %s",
        linked_count,
        len(feature_names),
        item_type,
        item_id,
    )


def find_or_create_feature_in_postgres(
    feature_name: str,
    feature_type: str,
    db: Optional[Session] = None,
) -> Optional[FeatureModel]:
    """
    Find or create a Feature record in Postgres.

    Public function for use by other modules.
    """
    if db is None:
        db_gen = get_db()
        session = next(db_gen)
        try:
            return _find_or_create_feature_impl(feature_name, feature_type, session)
        finally:
            session.close()
            next(db_gen, None)

    return _find_or_create_feature_impl(feature_name, feature_type, db)


def _find_or_create_feature_impl(
    feature_name: str,
    feature_type: str,
    db: Session,
) -> Optional[FeatureModel]:
    """Find or create a Feature record in Postgres."""
    # Normalize feature name (lowercase, strip whitespace)
    normalized_name = feature_name.strip().lower()

    # Try to find existing feature
    feature = (
        db.query(FeatureModel)
        .filter(FeatureModel.name.ilike(normalized_name))
        .filter(FeatureModel.feature_type == feature_type)
        .first()
    )

    if feature:
        return feature

    # Create new feature
    try:
        feature = FeatureModel(
            name=normalized_name,
            feature_type=feature_type,
        )
        db.add(feature)
        db.commit()
        db.refresh(feature)
        logger.debug(
            "[INGEST][FEATURES] Created new feature %s (ID: %s)",
            normalized_name,
            feature.id,
        )
        return feature
    except Exception as e:
        logger.error(
            "[INGEST][FEATURES] Error creating feature %s: %r",
            normalized_name,
            e,
        )
        db.rollback()
        return None


def _update_semantic_metadata_with_features(
    item_id: UUID,
    item_type: str,
    feature_mentions: List[str],
    db: Session,
) -> None:
    """Update semantic_metadata JSON field with feature mentions."""
    try:
        item: Optional[Union[Literature, Email]]
        if item_type == "literature":
            item = db.query(Literature).filter(Literature.id == item_id).first()
        elif item_type == "email":
            item = db.query(Email).filter(Email.id == item_id).first()
        else:
            logger.warning(
                "[INGEST][FEATURES] Cannot update semantic_metadata for item_type %s",
                item_type,
            )
            return

        if not item:
            logger.warning(
                "[INGEST][FEATURES] %s record %s not found",
                item_type,
                item_id,
            )
            return

        # Get existing metadata or create new
        semantic_meta = item.semantic_metadata or {}

        # Get existing feature mentions
        existing_mentions = semantic_meta.get("feature_mentions", [])

        # Merge with new mentions (avoid duplicates)
        all_mentions = sorted(set(existing_mentions + feature_mentions))

        # Update metadata
        semantic_meta["feature_mentions"] = all_mentions
        item.semantic_metadata = semantic_meta

        db.commit()
        logger.debug(
            "[INGEST][FEATURES] Updated semantic_metadata for %s %s with %d feature mentions",
            item_type,
            item_id,
            len(all_mentions),
        )
    except Exception as e:
        logger.error(
            "[INGEST][FEATURES] Error updating semantic_metadata for %s %s: %r",
            item_type,
            item_id,
            e,
        )
        db.rollback()


def batch_link_features_to_dataset_in_postgres(
    features: List[Tuple[str, str]],
    dataset_id: UUID,
    db: Optional[Session] = None,
    max_workers: int = 5,
) -> Dict[str, bool]:
    """
    Batch link features to a Postgres dataset with parallel processing.

    This function efficiently links multiple features to a dataset by:
    1. Processing features in parallel batches
    2. Finding or creating Feature records
    3. Linking features to the dataset via the many-to-many relationship

    Args:
        features: List of (feature_name, feature_type) tuples
        dataset_id: UUID of the target dataset
        db: Optional database session (creates new if None)
        max_workers: Maximum number of parallel workers (default: 5)

    Returns:
        Dictionary mapping feature_name -> success (True/False)

    Example:
        >>> features = [("PC(16:0/18:1)", "lipid"), ("PE(18:0/20:4)", "lipid")]
        >>> results = batch_link_features_to_dataset_in_postgres(
        ...     features=features,
        ...     dataset_id=dataset_uuid,
        ... )
        >>> results["PC(16:0/18:1)"]
        True
    """
    if db is None:
        db_gen = get_db()
        session = next(db_gen)
        try:
            return batch_link_features_to_dataset_in_postgres(
                features=features,
                dataset_id=dataset_id,
                db=session,
                max_workers=max_workers,
            )
        finally:
            session.close()
            next(db_gen, None)

    if not features:
        logger.debug(
            "[INGEST][FEATURES] No features to batch link for dataset %s",
            dataset_id,
        )
        return {}

    logger.info(
        "[INGEST][FEATURES] Batch linking %d feature(s) to dataset %s (max_workers=%d)",
        len(features),
        dataset_id,
        max_workers,
    )

    # Get dataset
    dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    if not dataset:
        logger.error(
            "[INGEST][FEATURES] Dataset %s not found",
            dataset_id,
        )
        return {name: False for name, _ in features}

    # OPTIMIZATION: Batch lookup all features first to avoid N+1 queries
    # Group features by type for efficient batch querying
    features_by_type: Dict[str, List[str]] = {}
    for name, ftype in features:
        ftype_key = ftype or ""
        if ftype_key not in features_by_type:
            features_by_type[ftype_key] = []
        features_by_type[ftype_key].append(name)

    # Batch lookup existing features (single query per type)
    feature_cache: Dict[Tuple[str, str], Optional[FeatureModel]] = {}
    existing_features: Dict[str, FeatureModel] = {}  # normalized_name -> FeatureModel

    for ftype, names in features_by_type.items():
        # Normalize all names
        normalized_names = [n.strip().lower() for n in names]

        # Single query to find all existing features of this type
        existing = (
            db.query(FeatureModel)
            .filter(FeatureModel.feature_type == ftype)
            .filter(FeatureModel.name.in_(normalized_names))
            .all()
        )

        # Build lookup map
        for feat in existing:
            if not feat.name:
                continue
            feat_type = feat.feature_type or ""
            existing_features[feat.name] = feat
            feature_cache[(feat.name, feat_type)] = feat

    # Create missing features in batch
    features_to_create: List[Tuple[str, str]] = []
    for name, ftype in features:
        normalized = name.strip().lower()
        ftype_key = ftype or ""
        if (normalized, ftype_key) not in feature_cache:
            features_to_create.append((normalized, ftype_key))

    if features_to_create:
        logger.debug(
            "[INGEST][FEATURES] Creating %d new feature(s) in batch",
            len(features_to_create),
        )
        for normalized, ftype in features_to_create:
            try:
                new_feature = FeatureModel(
                    name=normalized,
                    feature_type=ftype,
                )
                db.add(new_feature)
                feature_cache[(normalized, ftype)] = new_feature
            except Exception as e:
                logger.warning(
                    "[INGEST][FEATURES] Error creating feature %s (type: %s): %r",
                    normalized,
                    ftype,
                    e,
                )
                feature_cache[(normalized, ftype)] = None

        # Commit new features
        try:
            db.commit()
            # Refresh to get IDs
            for normalized, ftype in features_to_create:
                if feature_cache.get((normalized, ftype)):
                    db.refresh(feature_cache[(normalized, ftype)])
        except Exception as e:
            logger.error(
                "[INGEST][FEATURES] Error committing new features: %r",
                e,
            )
            db.rollback()
            # Mark failed features
            for normalized, ftype in features_to_create:
                feature_cache[(normalized, ftype)] = None

    # Results dictionary
    results: Dict[str, bool] = {}

    # Link features to dataset (no more database queries needed)
    for name, ftype in features:
        normalized = name.strip().lower()
        ftype_key = ftype or ""
        feature = feature_cache.get((normalized, ftype_key))

        if not feature:
            logger.warning(
                "[INGEST][FEATURES] Could not find/create feature %s (type: %s)",
                name,
                ftype,
            )
            results[name] = False
            continue

        # Link to dataset if not already linked
        if feature not in dataset.features:
            dataset.features.append(feature)
            logger.debug(
                "[INGEST][FEATURES] Linked feature %s to dataset %s",
                name,
                dataset_id,
            )
        else:
            logger.debug(
                "[INGEST][FEATURES] Feature %s already linked to dataset %s",
                name,
                dataset_id,
            )

        results[name] = True

    # Commit all changes
    try:
        db.commit()
        linked_count = sum(1 for v in results.values() if v)
        logger.info(
            "[INGEST][FEATURES] Batch linked %d/%d feature(s) to dataset %s",
            linked_count,
            len(features),
            dataset_id,
        )
    except Exception as e:
        logger.error(
            "[INGEST][FEATURES] Error committing feature links for dataset %s: %r",
            dataset_id,
            e,
        )
        db.rollback()
        # Mark all as failed
        results = {name: False for name, _ in features}

    return results
