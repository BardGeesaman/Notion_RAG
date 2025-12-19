"""
Dataset Comparison & Clustering.

Compares datasets across all omics types, computes similarity scores,
finds shared/differential features, and clusters datasets by similarity.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _get_dataset_name(dataset_id: str) -> str:
    """Get dataset name from Postgres."""
    with db_session() as db:
        try:
            dataset = None
            try:
                dataset_uuid = UUID(dataset_id)
                dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_uuid).first()
            except ValueError:
                dataset = db.query(DatasetModel).filter(DatasetModel.notion_page_id == dataset_id).first()

            if dataset and dataset.name:
                return dataset.name
            return f"Dataset {dataset_id[:8]}"
        except Exception as e:
            logger.debug("[ANALYSIS][DATASET-COMP] Could not fetch dataset name: %r", e)
            return f"Dataset {dataset_id[:8]}"


@dataclass
class DatasetComparison:
    """
    Results from comparing two datasets.

    Attributes:
        dataset1_id: Notion page ID of first dataset
        dataset2_id: Notion page ID of second dataset
        dataset1_name: Name of first dataset
        dataset2_name: Name of second dataset
        overall_similarity: Overall similarity score (0-1)
        similarity_by_omics: Similarity scores by omics type
        shared_features: Features present in both datasets (by omics type)
        unique_to_dataset1: Features only in dataset1 (by omics type)
        unique_to_dataset2: Features only in dataset2 (by omics type)
        jaccard_similarity: Jaccard similarity coefficient
    """
    dataset1_id: str
    dataset2_id: str
    dataset1_name: str
    dataset2_name: str
    overall_similarity: float
    similarity_by_omics: Dict[str, float] = field(default_factory=dict)
    shared_features: Dict[str, Set[str]] = field(default_factory=dict)
    unique_to_dataset1: Dict[str, Set[str]] = field(default_factory=dict)
    unique_to_dataset2: Dict[str, Set[str]] = field(default_factory=dict)
    jaccard_similarity: float = 0.0


@dataclass
class DatasetCluster:
    """
    Represents a cluster of similar datasets.

    Attributes:
        cluster_id: Unique cluster identifier
        dataset_ids: List of dataset page IDs in this cluster
        dataset_names: List of dataset names
        average_similarity: Average pairwise similarity within cluster
        representative_dataset: Most representative dataset ID
    """
    cluster_id: int
    dataset_ids: List[str] = field(default_factory=list)
    dataset_names: List[str] = field(default_factory=list)
    average_similarity: float = 0.0
    representative_dataset: Optional[str] = None


def compute_jaccard_similarity(set1: Set[str], set2: Set[str]) -> float:
    """
    Compute Jaccard similarity coefficient between two sets.

    Jaccard = |A ∩ B| / |A ∪ B|

    Args:
        set1: First set
        set2: Second set

    Returns:
        Jaccard similarity (0-1)
    """
    if not set1 and not set2:
        return 1.0  # Both empty = identical

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    if union == 0:
        return 0.0

    return intersection / union


def compare_datasets(
    dataset1_id: str,
    dataset2_id: str,
    use_cache: bool = True,
) -> DatasetComparison:
    """
    Compare two datasets across all omics types.

    Args:
        dataset1_id: Notion page ID of first dataset (with dashes)
        dataset2_id: Notion page ID of second dataset (with dashes)
        use_cache: Whether to use feature cache

    Returns:
        DatasetComparison object
    """
    logger.info(
        "[ANALYSIS][DATASET-COMP] Comparing datasets %s and %s",
        dataset1_id,
        dataset2_id,
    )

    # Get dataset names from Postgres
    dataset1_name = _get_dataset_name(dataset1_id)
    dataset2_name = _get_dataset_name(dataset2_id)

    # Extract features from both datasets
    features1 = extract_dataset_features_by_type(dataset1_id, use_cache=use_cache)
    features2 = extract_dataset_features_by_type(dataset2_id, use_cache=use_cache)

    # Compare by omics type
    similarity_by_omics: Dict[str, float] = {}
    shared_features: Dict[str, Set[str]] = {}
    unique_to_dataset1: Dict[str, Set[str]] = {}
    unique_to_dataset2: Dict[str, Set[str]] = {}

    all_omics_types = set(features1.keys()) | set(features2.keys())

    for omics_type in all_omics_types:
        set1 = features1.get(omics_type, set())
        set2 = features2.get(omics_type, set())

        # Jaccard similarity
        similarity = compute_jaccard_similarity(set1, set2)
        similarity_by_omics[omics_type] = similarity

        # Shared and unique features
        shared = set1 & set2
        unique1 = set1 - set2
        unique2 = set2 - set1

        if shared:
            shared_features[omics_type] = shared
        if unique1:
            unique_to_dataset1[omics_type] = unique1
        if unique2:
            unique_to_dataset2[omics_type] = unique2

    # Overall similarity (weighted average by feature count)
    total_features1 = sum(len(f) for f in features1.values())
    total_features2 = sum(len(f) for f in features2.values())
    total_features = total_features1 + total_features2

    if total_features > 0:
        overall_similarity = sum(
            similarity_by_omics.get(omics_type, 0.0) * (len(features1.get(omics_type, set())) + len(features2.get(omics_type, set())))
            for omics_type in all_omics_types
        ) / total_features
    else:
        overall_similarity = 0.0

    # Overall Jaccard similarity (all features combined)
    all_features1 = set()
    all_features2 = set()
    for features in features1.values():
        all_features1.update(features)
    for features in features2.values():
        all_features2.update(features)

    jaccard_similarity = compute_jaccard_similarity(all_features1, all_features2)

    return DatasetComparison(
        dataset1_id=dataset1_id,
        dataset2_id=dataset2_id,
        dataset1_name=dataset1_name,
        dataset2_name=dataset2_name,
        overall_similarity=overall_similarity,
        similarity_by_omics=similarity_by_omics,
        shared_features=shared_features,
        unique_to_dataset1=unique_to_dataset1,
        unique_to_dataset2=unique_to_dataset2,
        jaccard_similarity=jaccard_similarity,
    )


def compare_multiple_datasets(
    dataset_ids: List[str],
    use_cache: bool = True,
) -> List[DatasetComparison]:
    """
    Compare multiple datasets pairwise.

    Args:
        dataset_ids: List of dataset page IDs
        use_cache: Whether to use feature cache

    Returns:
        List of DatasetComparison objects (all pairwise comparisons)
    """
    logger.info(
        "[ANALYSIS][DATASET-COMP] Comparing %d datasets pairwise",
        len(dataset_ids),
    )

    comparisons: List[DatasetComparison] = []

    # Compare all pairs
    for i in range(len(dataset_ids)):
        for j in range(i + 1, len(dataset_ids)):
            try:
                comparison = compare_datasets(
                    dataset_ids[i],
                    dataset_ids[j],
                    use_cache=use_cache,
                )
                comparisons.append(comparison)
            except Exception as e:
                logger.warning(
                    "[ANALYSIS][DATASET-COMP] Error comparing %s and %s: %r",
                    dataset_ids[i],
                    dataset_ids[j],
                    e,
                )
                continue

    logger.info(
        "[ANALYSIS][DATASET-COMP] Completed %d pairwise comparisons",
        len(comparisons),
    )

    return comparisons


def cluster_datasets_by_similarity(
    comparisons: List[DatasetComparison],
    similarity_threshold: float = 0.5,
) -> List[DatasetCluster]:
    """
    Cluster datasets by similarity using a simple threshold-based approach.

    Args:
        comparisons: List of DatasetComparison objects
        similarity_threshold: Minimum similarity to be in same cluster

    Returns:
        List of DatasetCluster objects
    """
    logger.info(
        "[ANALYSIS][DATASET-COMP] Clustering datasets (threshold=%.2f)",
        similarity_threshold,
    )

    # Build similarity graph
    dataset_similarities: Dict[str, Dict[str, float]] = defaultdict(dict)

    for comp in comparisons:
        dataset_similarities[comp.dataset1_id][comp.dataset2_id] = comp.overall_similarity
        dataset_similarities[comp.dataset2_id][comp.dataset1_id] = comp.overall_similarity

    # Simple clustering: connected components where similarity >= threshold
    all_dataset_ids = set(dataset_similarities.keys())
    visited: Set[str] = set()
    clusters: List[DatasetCluster] = []
    cluster_id = 0

    for dataset_id in all_dataset_ids:
        if dataset_id in visited:
            continue

        # Start new cluster
        cluster_dataset_ids = [dataset_id]
        visited.add(dataset_id)

        # Find all connected datasets (similarity >= threshold)
        queue = [dataset_id]
        while queue:
            current_id = queue.pop(0)
            for neighbor_id, similarity in dataset_similarities[current_id].items():
                if neighbor_id not in visited and similarity >= similarity_threshold:
                    cluster_dataset_ids.append(neighbor_id)
                    visited.add(neighbor_id)
                    queue.append(neighbor_id)

        # Calculate average similarity within cluster
        if len(cluster_dataset_ids) > 1:
            similarities = []
            for i in range(len(cluster_dataset_ids)):
                for j in range(i + 1, len(cluster_dataset_ids)):
                    id1 = cluster_dataset_ids[i]
                    id2 = cluster_dataset_ids[j]
                    if id2 in dataset_similarities[id1]:
                        similarities.append(dataset_similarities[id1][id2])

            avg_similarity = sum(similarities) / len(similarities) if similarities else 0.0
        else:
            avg_similarity = 1.0

        # Get dataset names
        dataset_names = []
        for d_id in cluster_dataset_ids:
            # Try to get name from comparisons
            name = f"Dataset {d_id[:8]}"
            for comp in comparisons:
                if comp.dataset1_id == d_id:
                    name = comp.dataset1_name
                    break
                elif comp.dataset2_id == d_id:
                    name = comp.dataset2_name
                    break
            dataset_names.append(name)

        # Representative dataset (one with highest average similarity to others)
        if len(cluster_dataset_ids) > 1:
            representative = max(
                cluster_dataset_ids,
                key=lambda d_id: sum(
                    dataset_similarities[d_id].get(other_id, 0.0)
                    for other_id in cluster_dataset_ids
                    if other_id != d_id
                ),
            )
        else:
            representative = cluster_dataset_ids[0]

        clusters.append(
            DatasetCluster(
                cluster_id=cluster_id,
                dataset_ids=cluster_dataset_ids,
                dataset_names=dataset_names,
                average_similarity=avg_similarity,
                representative_dataset=representative,
            )
        )
        cluster_id += 1

    logger.info(
        "[ANALYSIS][DATASET-COMP] Found %d clusters from %d datasets",
        len(clusters),
        len(all_dataset_ids),
    )

    return clusters


def generate_comparison_report(
    comparison: DatasetComparison,
) -> str:
    """
    Generate a human-readable comparison report.

    Args:
        comparison: DatasetComparison object

    Returns:
        Markdown-formatted report
    """
    report = "# Dataset Comparison Report\n\n"
    report += f"## {comparison.dataset1_name} vs {comparison.dataset2_name}\n\n"

    report += "### Overall Similarity\n\n"
    report += f"- **Overall Similarity**: {comparison.overall_similarity:.3f}\n"
    report += f"- **Jaccard Similarity**: {comparison.jaccard_similarity:.3f}\n\n"

    if comparison.similarity_by_omics:
        report += "### Similarity by Omics Type\n\n"
        for omics_type, similarity in sorted(comparison.similarity_by_omics.items()):
            report += f"- **{omics_type.capitalize()}**: {similarity:.3f}\n"
        report += "\n"

    if comparison.shared_features:
        report += "### Shared Features\n\n"
        total_shared = sum(len(f) for f in comparison.shared_features.values())
        report += f"**Total shared features**: {total_shared}\n\n"

        for omics_type, features in sorted(comparison.shared_features.items()):
            report += f"#### {omics_type.capitalize()} ({len(features)} features)\n\n"
            # Show first 20 features
            feature_list = sorted(list(features))[:20]
            report += f"{', '.join(feature_list)}"
            if len(features) > 20:
                report += f" (and {len(features) - 20} more)"
            report += "\n\n"

    if comparison.unique_to_dataset1:
        report += f"### Unique to {comparison.dataset1_name}\n\n"
        total_unique1 = sum(len(f) for f in comparison.unique_to_dataset1.values())
        report += f"**Total unique features**: {total_unique1}\n\n"

        for omics_type, features in sorted(comparison.unique_to_dataset1.items()):
            report += f"#### {omics_type.capitalize()} ({len(features)} features)\n\n"
            feature_list = sorted(list(features))[:20]
            report += f"{', '.join(feature_list)}"
            if len(features) > 20:
                report += f" (and {len(features) - 20} more)"
            report += "\n\n"

    if comparison.unique_to_dataset2:
        report += f"### Unique to {comparison.dataset2_name}\n\n"
        total_unique2 = sum(len(f) for f in comparison.unique_to_dataset2.values())
        report += f"**Total unique features**: {total_unique2}\n\n"

        for omics_type, features in sorted(comparison.unique_to_dataset2.items()):
            report += f"#### {omics_type.capitalize()} ({len(features)} features)\n\n"
            feature_list = sorted(list(features))[:20]
            report += f"{', '.join(feature_list)}"
            if len(features) > 20:
                report += f" (and {len(features) - 20} more)"
            report += "\n\n"

    return report


def generate_clustering_report(
    clusters: List[DatasetCluster],
) -> str:
    """
    Generate a human-readable clustering report.

    Args:
        clusters: List of DatasetCluster objects

    Returns:
        Markdown-formatted report
    """
    report = "# Dataset Clustering Report\n\n"
    report += f"Found **{len(clusters)}** cluster(s) from dataset comparisons.\n\n"

    for cluster in clusters:
        report += f"## Cluster {cluster.cluster_id + 1}\n\n"
        report += f"- **Number of datasets**: {len(cluster.dataset_ids)}\n"
        report += f"- **Average similarity**: {cluster.average_similarity:.3f}\n"
        report += f"- **Representative dataset**: {cluster.representative_dataset}\n\n"

        report += "### Datasets in Cluster\n\n"
        for i, (dataset_id, dataset_name) in enumerate(zip(cluster.dataset_ids, cluster.dataset_names), 1):
            report += f"{i}. {dataset_name} (`{dataset_id}`)\n"
        report += "\n"

    return report
