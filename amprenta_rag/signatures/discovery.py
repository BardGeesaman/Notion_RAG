"""
Automated signature discovery from datasets.

Analyzes multiple datasets to discover recurring patterns and generate
candidate signatures based on feature co-occurrence, direction consistency,
and statistical significance.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class FeatureOccurrence:
    """
    Represents a feature occurrence in a dataset.

    Attributes:
        feature_name: Normalized feature name
        feature_type: Omics type (gene, protein, metabolite, lipid)
        direction: Direction of change (↑, ↓, or None)
        dataset_id: Dataset page ID where this feature was found
        fold_change: Optional fold change value
        p_value: Optional p-value
    """
    feature_name: str
    feature_type: str
    direction: Optional[str] = None
    dataset_id: str = ""
    fold_change: Optional[float] = None
    p_value: Optional[float] = None


@dataclass
class SignatureCandidate:
    """
    Represents a candidate signature discovered from datasets.

    Attributes:
        name: Suggested signature name
        features: List of features in the signature
        feature_types: Set of omics types included
        datasets: List of dataset IDs where this pattern was found
        co_occurrence_score: How often features co-occur
        direction_consistency: Consistency of direction across datasets
        support_count: Number of datasets supporting this signature
        confidence: Overall confidence score (0-1)
    """
    name: str
    features: List[FeatureOccurrence] = field(default_factory=list)
    feature_types: Set[str] = field(default_factory=set)
    datasets: List[str] = field(default_factory=list)
    co_occurrence_score: float = 0.0
    direction_consistency: float = 0.0
    support_count: int = 0
    confidence: float = 0.0


def extract_features_from_datasets(
    dataset_ids: List[str],
    min_features: int = 3,
    use_parallel: bool = False,
    max_workers: int = 5,
) -> Dict[str, List[FeatureOccurrence]]:
    """
    Extract features from multiple datasets.

    Args:
        dataset_ids: List of dataset page IDs
        min_features: Minimum number of features per dataset to include
        use_parallel: Whether to extract features in parallel (default: False)
        max_workers: Maximum number of parallel workers (if use_parallel=True)

    Returns:
        Dictionary mapping dataset_id to list of FeatureOccurrence
    """
    from amprenta_rag.ingestion.multi_omics_scoring import (
        extract_dataset_features_by_type,
    )

    logger.info(
        "[SIG-DISC] Extracting features from %d datasets%s",
        len(dataset_ids),
        " (parallel)" if use_parallel else "",
    )

    all_features: Dict[str, List[FeatureOccurrence]] = {}

    def extract_one_dataset(dataset_id: str) -> tuple[str, List[FeatureOccurrence]]:
        """Extract features from a single dataset."""
        try:
            features_by_type = extract_dataset_features_by_type(
                dataset_id,
                use_cache=True,
            )

            occurrences: List[FeatureOccurrence] = []
            for omics_type, feature_set in features_by_type.items():
                if not feature_set:
                    continue
                for feature_name in feature_set:
                    occurrence = FeatureOccurrence(
                        feature_name=feature_name,
                        feature_type=omics_type,
                        dataset_id=dataset_id,
                    )
                    occurrences.append(occurrence)

            return (dataset_id, occurrences)
        except Exception as e:
            logger.warning(
                "[SIG-DISC] Error extracting features from dataset %s: %r",
                dataset_id,
                e,
            )
            return (dataset_id, [])

    if use_parallel:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from tqdm import tqdm

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(extract_one_dataset, dataset_id): dataset_id
                for dataset_id in dataset_ids
            }

            # Use tqdm if available for progress
            if tqdm:
                futures_iter = tqdm(
                    as_completed(futures),
                    total=len(dataset_ids),
                    desc="Extracting features",
                    unit="dataset",
                )
            else:
                futures_iter = as_completed(futures)

            for future in futures_iter:
                dataset_id, occurrences = future.result()
                if len(occurrences) >= min_features:
                    all_features[dataset_id] = occurrences
    else:
        for dataset_id in dataset_ids:
            _, occurrences = extract_one_dataset(dataset_id)
            if len(occurrences) >= min_features:
                all_features[dataset_id] = occurrences
            elif len(occurrences) > 0:
                logger.debug(
                    "[SIG-DISC] Dataset %s has only %d features (min: %d), skipping",
                    dataset_id,
                    len(occurrences),
                    min_features,
                )

    total_features = sum(len(feats) for feats in all_features.values())
    logger.info(
        "[SIG-DISC] Extracted %d total features from %d datasets",
        total_features,
        len(all_features),
    )

    return all_features


def compute_feature_co_occurrence(
    all_features: Dict[str, List[FeatureOccurrence]],
    min_co_occurrence: int = 2,
) -> Dict[Tuple[str, str], int]:
    """
    Compute co-occurrence counts for feature pairs.

    Args:
        all_features: Dictionary mapping dataset_id to features
        min_co_occurrence: Minimum co-occurrence count to include

    Returns:
        Dictionary mapping (feature1, feature2) tuples to co-occurrence count
    """
    logger.info("[SIG-DISC] Computing feature co-occurrence...")

    # Build feature sets per dataset
    dataset_feature_sets: Dict[str, Set[str]] = {}
    for dataset_id, features in all_features.items():
        feature_set = {f.feature_name for f in features}
        dataset_feature_sets[dataset_id] = feature_set

    # Count co-occurrences
    co_occurrence_counts: Dict[Tuple[str, str], int] = defaultdict(int)

    for dataset_id, feature_set in dataset_feature_sets.items():
        features_list = sorted(list(feature_set))
        for i, feat1 in enumerate(features_list):
            for feat2 in features_list[i + 1:]:
                # Use canonical ordering (alphabetical) to avoid duplicates
                pair = (feat1, feat2) if feat1 < feat2 else (feat2, feat1)
                co_occurrence_counts[pair] += 1

    # Filter by minimum co-occurrence
    filtered_counts = {
        pair: count
        for pair, count in co_occurrence_counts.items()
        if count >= min_co_occurrence
    }

    logger.info(
        "[SIG-DISC] Found %d feature pairs with co-occurrence >= %d",
        len(filtered_counts),
        min_co_occurrence,
    )

    return filtered_counts


def compute_direction_consistency(
    feature_name: str,
    datasets: List[str],
    all_features: Dict[str, List[FeatureOccurrence]],
) -> float:
    """
    Compute direction consistency for a feature across datasets.

    Args:
        feature_name: Feature to check
        datasets: List of dataset IDs where feature appears
        all_features: Dictionary mapping dataset_id to features

    Returns:
        Consistency score (0-1), where 1.0 means all datasets agree on direction
    """
    directions: List[str] = []

    for dataset_id in datasets:
        features = all_features.get(dataset_id, [])
        for feat in features:
            if feat.feature_name == feature_name and feat.direction:
                directions.append(feat.direction)

    if not directions:
        return 0.0

    # Count direction occurrences
    direction_counts = Counter(directions)
    most_common_count = direction_counts.most_common(1)[0][1]

    # Consistency = fraction of datasets with the most common direction
    consistency = most_common_count / len(directions)

    return consistency


def discover_signature_candidates(
    dataset_ids: List[str],
    min_support: int = 3,
    min_features: int = 3,
    max_features: int = 20,
    min_co_occurrence: int = 2,
    min_confidence: float = 0.5,
) -> List[SignatureCandidate]:
    """
    Discover candidate signatures from datasets.

    Args:
        dataset_ids: List of dataset page IDs to analyze
        min_support: Minimum number of datasets that must contain the pattern
        min_features: Minimum number of features in a signature
        max_features: Maximum number of features in a signature
        min_co_occurrence: Minimum co-occurrence count for feature pairs
        min_confidence: Minimum confidence score to include candidate

    Returns:
        List of SignatureCandidate objects, sorted by confidence
    """
    logger.info(
        "[SIG-DISC] Discovering signature candidates from %d datasets",
        len(dataset_ids),
    )

    # Extract features from all datasets (sequential for now, can be parallelized)
    all_features = extract_features_from_datasets(
        dataset_ids,
        min_features=1,
        use_parallel=False,  # Can be made configurable
    )

    if not all_features:
        logger.warning("[SIG-DISC] No features extracted from datasets")
        return []

    # Compute co-occurrence
    co_occurrence = compute_feature_co_occurrence(
        all_features,
        min_co_occurrence=min_co_occurrence,
    )

    if not co_occurrence:
        logger.warning("[SIG-DISC] No significant co-occurrences found")
        return []

    # Build feature clusters based on co-occurrence
    # Use a simple clustering approach: features that co-occur frequently
    # form clusters

    # Build feature graph (adjacency based on co-occurrence)
    feature_graph: Dict[str, Set[str]] = defaultdict(set)
    feature_support: Dict[str, Set[str]] = defaultdict(set)  # datasets containing feature

    for (feat1, feat2), count in co_occurrence.items():
        feature_graph[feat1].add(feat2)
        feature_graph[feat2].add(feat1)

        # Track which datasets contain these features
        for dataset_id, features in all_features.items():
            feature_names = {f.feature_name for f in features}
            if feat1 in feature_names:
                feature_support[feat1].add(dataset_id)
            if feat2 in feature_names:
                feature_support[feat2].add(dataset_id)

    # Find clusters (connected components in the graph)
    visited: Set[str] = set()
    clusters: List[Set[str]] = []

    def dfs_cluster(feature: str, cluster: Set[str]) -> None:
        """Depth-first search to find connected features."""
        if feature in visited:
            return
        visited.add(feature)
        cluster.add(feature)

        for neighbor in feature_graph.get(feature, set()):
            if neighbor not in visited:
                dfs_cluster(neighbor, cluster)

    for feature in feature_graph:
        if feature not in visited:
            cluster: Set[str] = set()
            dfs_cluster(feature, cluster)
            if len(cluster) >= min_features:
                clusters.append(cluster)

    logger.info(
        "[SIG-DISC] Found %d feature clusters (size >= %d)",
        len(clusters),
        min_features,
    )

    # Convert clusters to signature candidates
    candidates: List[SignatureCandidate] = []

    for i, cluster in enumerate(clusters):
        if len(cluster) > max_features:
            # Take top features by support count
            cluster_list = sorted(
                cluster,
                key=lambda f: len(feature_support.get(f, set())),
                reverse=True,
            )
            cluster = set(cluster_list[:max_features])

        # Find datasets that contain all features in cluster
        cluster_features = list(cluster)
        if not cluster_features:
            continue

        # Find intersection of datasets containing all features
        supporting_datasets = feature_support.get(cluster_features[0], set())
        for feat in cluster_features[1:]:
            supporting_datasets &= feature_support.get(feat, set())

        if len(supporting_datasets) < min_support:
            continue

        # Build feature occurrences
        feature_occurrences: List[FeatureOccurrence] = []
        feature_types: Set[str] = set()

        for feat_name in cluster_features:
            # Find feature type from first occurrence
            feat_type = None
            for dataset_id, features in all_features.items():
                for feat in features:
                    if feat.feature_name == feat_name:
                        feat_type = feat.feature_type
                        break
                if feat_type:
                    break

            if not feat_type:
                continue

            # Compute direction consistency
            direction_consistency = compute_direction_consistency(
                feat_name,
                list(supporting_datasets),
                all_features,
            )

            # Determine most common direction
            directions: List[str] = []
            for dataset_id in supporting_datasets:
                features = all_features.get(dataset_id, [])
                for feat in features:
                    if feat.feature_name == feat_name and feat.direction:
                        directions.append(feat.direction)

            most_common_direction = None
            if directions:
                direction_counts = Counter(directions)
                most_common_direction = direction_counts.most_common(1)[0][0]

            occurrence = FeatureOccurrence(
                feature_name=feat_name,
                feature_type=feat_type,
                direction=most_common_direction,
                dataset_id="",  # Will be set per dataset
            )
            feature_occurrences.append(occurrence)
            feature_types.add(feat_type)

        if not feature_occurrences:
            continue

        # Compute scores
        support_count = len(supporting_datasets)
        co_occurrence_score = len(supporting_datasets) / len(dataset_ids)

        # Average direction consistency
        direction_consistencies = [
            compute_direction_consistency(
                feat.feature_name,
                list(supporting_datasets),
                all_features,
            )
            for feat in feature_occurrences
        ]
        avg_direction_consistency = (
            sum(direction_consistencies) / len(direction_consistencies)
            if direction_consistencies
            else 0.0
        )

        # Confidence = weighted combination of support and consistency
        confidence = (
            0.6 * co_occurrence_score + 0.4 * avg_direction_consistency
        )

        if confidence < min_confidence:
            continue

        # Generate signature name
        signature_name = f"Discovered-Signature-{i+1}"
        if feature_types:
            types_str = "-".join(sorted(feature_types))
            signature_name = f"Discovered-{types_str}-{i+1}"

        candidate = SignatureCandidate(
            name=signature_name,
            features=feature_occurrences,
            feature_types=feature_types,
            datasets=list(supporting_datasets),
            co_occurrence_score=co_occurrence_score,
            direction_consistency=avg_direction_consistency,
            support_count=support_count,
            confidence=confidence,
        )
        candidates.append(candidate)

    # Sort by confidence
    candidates.sort(key=lambda c: c.confidence, reverse=True)

    logger.info(
        "[SIG-DISC] Discovered %d signature candidates (confidence >= %.2f)",
        len(candidates),
        min_confidence,
    )

    return candidates


def export_candidates_to_tsv(
    candidates: List[SignatureCandidate],
    output_path: str,
) -> None:
    """
    Export signature candidates to TSV file.

    Args:
        candidates: List of SignatureCandidate objects
        output_path: Path to output TSV file
    """
    from pathlib import Path

    logger.info(
        "[SIG-DISC] Exporting %d candidates to %s",
        len(candidates),
        output_path,
    )

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        # Write header
        f.write("signature_name\tfeature_type\tfeature_name\tdirection\tweight\tsupport_count\tconfidence\n")

        # Write candidates
        for candidate in candidates:
            # Assign weights based on confidence and support
            base_weight = candidate.confidence

            for feat in candidate.features:
                # Weight can be adjusted based on feature importance
                weight = base_weight

                f.write(
                    f"{candidate.name}\t"
                    f"{feat.feature_type}\t"
                    f"{feat.feature_name}\t"
                    f"{feat.direction or ''}\t"
                    f"{weight:.3f}\t"
                    f"{candidate.support_count}\t"
                    f"{candidate.confidence:.3f}\n"
                )

    logger.info("[SIG-DISC] Exported candidates to %s", output_path)


def export_candidates_to_json(
    candidates: List[SignatureCandidate],
    output_path: str,
) -> None:
    """
    Export signature candidates to JSON file.

    Args:
        candidates: List of SignatureCandidate objects
        output_path: Path to output JSON file
    """
    import json
    from pathlib import Path

    logger.info(
        "[SIG-DISC] Exporting %d candidates to %s (JSON)",
        len(candidates),
        output_path,
    )

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Convert candidates to JSON-serializable format
    export_data = {
        "candidates": [
            {
                "name": c.name,
                "features": [
                    {
                        "feature_name": f.feature_name,
                        "feature_type": f.feature_type,
                        "direction": f.direction,
                    }
                    for f in c.features
                ],
                "feature_types": list(c.feature_types),
                "datasets": c.datasets,
                "co_occurrence_score": c.co_occurrence_score,
                "direction_consistency": c.direction_consistency,
                "support_count": c.support_count,
                "confidence": c.confidence,
            }
            for c in candidates
        ],
        "summary": {
            "total_candidates": len(candidates),
            "total_features": sum(len(c.features) for c in candidates),
            "avg_confidence": (
                sum(c.confidence for c in candidates) / len(candidates)
                if candidates
                else 0.0
            ),
        },
    }

    with open(output_file, "w") as f:
        json.dump(export_data, f, indent=2)

    logger.info("[SIG-DISC] Exported candidates to %s", output_path)

