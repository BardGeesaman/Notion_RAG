"""Pattern detection utilities for cross-dataset analysis."""
from __future__ import annotations

from typing import List, Dict, Set, Any
from uuid import UUID

from amprenta_rag.database.models import Dataset, Signature, SignatureComponent, Feature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def find_recurring_features(
    dataset_ids: List[str],
    db,
    min_occurrence: int = 2,
) -> List[Dict[str, any]]:
    """
    Find features that appear across multiple datasets.

    Args:
        dataset_ids: List of dataset UUIDs
        db: Database session
        min_occurrence: Minimum number of datasets a feature must appear in

    Returns:
        List of dicts with feature_name, occurrence_count, datasets
    """
    if not dataset_ids:
        logger.warning("[PATTERN] Empty dataset_ids list")
        return []

    try:
        # Get signatures linked to these datasets
        dataset_uuids = [UUID(did) for did in dataset_ids]
        signatures = (
            db.query(Signature)
            .join(Signature.datasets)
            .filter(Dataset.id.in_(dataset_uuids))
            .distinct()
            .all()
        )

        logger.info("[PATTERN] Found %d signatures across %d datasets", len(signatures), len(dataset_ids))

        # Collect features from signature components
        feature_datasets: Dict[str, Set[str]] = {}  # feature_name -> set of dataset_ids

        for signature in signatures:
            # Get datasets linked to this signature
            sig_datasets = {str(d.id) for d in signature.datasets if str(d.id) in dataset_ids}

            # Get components
            components = db.query(SignatureComponent).filter(
                SignatureComponent.signature_id == signature.id
            ).all()

            for component in components:
                # Get feature name (from component.feature_name or linked Feature)
                if component.feature_name:
                    feature_name = component.feature_name
                elif component.feature_id:
                    feature = db.query(Feature).filter(Feature.id == component.feature_id).first()
                    feature_name = feature.name if feature else None
                else:
                    continue

                if feature_name:
                    if feature_name not in feature_datasets:
                        feature_datasets[feature_name] = set()
                    feature_datasets[feature_name].update(sig_datasets)

        # Filter by min_occurrence and format results
        results = []
        for feature_name, dataset_set in feature_datasets.items():
            occurrence_count = len(dataset_set)
            if occurrence_count >= min_occurrence:
                results.append({
                    "feature_name": feature_name,
                    "occurrence_count": occurrence_count,
                    "datasets": sorted(list(dataset_set)),
                })

        # Sort by occurrence count descending
        results.sort(key=lambda x: x["occurrence_count"], reverse=True)

        logger.info("[PATTERN] Found %d recurring features (min_occurrence=%d)", len(results), min_occurrence)
        return results

    except Exception as e:
        logger.error("[PATTERN] Error finding recurring features: %r", e)
        return []


def calculate_overlap(dataset_ids: List[str], db) -> Dict[str, Set[str]]:
    """
    Calculate feature overlap between datasets.

    Args:
        dataset_ids: List of dataset UUIDs
        db: Database session

    Returns:
        Dict mapping dataset_id to set of feature names
    """
    if not dataset_ids:
        logger.warning("[PATTERN] Empty dataset_ids list")
        return {}

    try:
        [UUID(did) for did in dataset_ids]
        overlap: Dict[str, Set[str]] = {}

        for dataset_id in dataset_ids:
            dataset_uuid = UUID(dataset_id)
            dataset = db.query(Dataset).filter(Dataset.id == dataset_uuid).first()

            if not dataset:
                continue

            feature_set: Set[str] = set()

            # Get features directly linked to dataset
            for feature in dataset.features:
                feature_set.add(feature.name)

            # Get features from signatures linked to this dataset
            for signature in dataset.signatures:
                components = db.query(SignatureComponent).filter(
                    SignatureComponent.signature_id == signature.id
                ).all()

                for component in components:
                    if component.feature_name:
                        feature_set.add(component.feature_name)
                    elif component.feature_id:
                        feature = db.query(Feature).filter(Feature.id == component.feature_id).first()
                        if feature:
                            feature_set.add(feature.name)

            overlap[dataset_id] = feature_set

        logger.info("[PATTERN] Calculated overlap for %d datasets", len(overlap))
        return overlap

    except Exception as e:
        logger.error("[PATTERN] Error calculating overlap: %r", e)
        return {}


def get_cross_dataset_summary(dataset_ids: List[str], db) -> Dict[str, Any]:
    """
    Get summary statistics for features across datasets.

    Args:
        dataset_ids: List of dataset UUIDs
        db: Database session

    Returns:
        Dict with:
        - total_features_per_dataset: {dataset_id: count}
        - overlap_counts: {dataset_pair: overlap_count}
        - top_recurring_features: List of recurring features
        - unique_features_per_dataset: {dataset_id: count}
    """
    if not dataset_ids:
        logger.warning("[PATTERN] Empty dataset_ids list")
        return {
            "total_features_per_dataset": {},
            "overlap_counts": {},
            "top_recurring_features": [],
            "unique_features_per_dataset": {},
        }

    try:
        # Get feature overlap
        overlap = calculate_overlap(dataset_ids, db)

        # Calculate per-dataset stats
        total_features_per_dataset = {
            did: len(features) for did, features in overlap.items()
        }

        # Calculate pairwise overlaps
        overlap_counts: Dict[str, int] = {}
        dataset_list = list(overlap.keys())
        for i, did1 in enumerate(dataset_list):
            for did2 in dataset_list[i+1:]:
                overlap_size = len(overlap[did1] & overlap[did2])
                if overlap_size > 0:
                    pair_key = f"{did1[:8]}...-{did2[:8]}..."
                    overlap_counts[pair_key] = overlap_size

        # Get recurring features
        recurring = find_recurring_features(dataset_ids, db, min_occurrence=2)
        top_recurring = recurring[:10]  # Top 10

        # Calculate unique features (features only in this dataset)
        unique_features_per_dataset = {}
        for did, features in overlap.items():
            other_features = set()
            for other_did, other_feat_set in overlap.items():
                if other_did != did:
                    other_features.update(other_feat_set)
            unique = features - other_features
            unique_features_per_dataset[did] = len(unique)

        summary = {
            "total_features_per_dataset": total_features_per_dataset,
            "overlap_counts": overlap_counts,
            "top_recurring_features": top_recurring,
            "unique_features_per_dataset": unique_features_per_dataset,
        }

        logger.info("[PATTERN] Generated cross-dataset summary for %d datasets", len(dataset_ids))
        return summary

    except Exception as e:
        logger.error("[PATTERN] Error generating summary: %r", e)
        return {
            "total_features_per_dataset": {},
            "overlap_counts": {},
            "top_recurring_features": [],
            "unique_features_per_dataset": {},
        }
