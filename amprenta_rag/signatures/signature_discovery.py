"""
Automated Signature Discovery Algorithm (v1 Experimental)

This module implements the core algorithm for discovering candidate multi-omics signatures
by analyzing feature co-occurrence patterns across multiple datasets.

**Algorithm Overview**:
1. Build feature-to-datasets mapping (inverted index)
2. Filter features by minimum support (must appear in N datasets)
3. Cluster features with similar co-occurrence patterns (overlap-based clustering)
4. Generate candidate signatures from clusters

**Version**: 1.0 (Experimental)
**Status**: For exploratory use - requires validation before production

**Limitations**:
- Does NOT extract or use direction information (all components have direction=None)
- Does NOT perform statistical significance testing
- Simple overlap-based clustering (not ML-based)
- Requires pre-normalized feature names
- No temporal or causal analysis

**See**: docs/SIGNATURE_DISCOVERY.md for complete documentation and validation workflow

**Future Enhancements** (v2.0+):
- Direction extraction from fold-change data
- Statistical significance testing (Fisher's exact, permutation)
- Cross-dataset normalization
- ML-based clustering algorithms
"""

import hashlib
from collections import defaultdict
from typing import List

from amprenta_rag.domain.signatures_discovery import Component, DiscoveredSignature, DiscoveryDatasetSummary


def discover_signatures_from_datasets(
    datasets: List[DiscoveryDatasetSummary],
    min_support: int = 2,
    min_overlap: float = 0.3,
) -> List[DiscoveredSignature]:
    """
    Discover candidate signatures from datasets using co-occurrence analysis.

    **Algorithm**:
    1. Build inverted index: feature → set of datasets containing it
    2. Filter features: keep only those appearing in ≥ min_support datasets
    3. Cluster features: group features with ≥ min_overlap co-occurrence ratio
    4. Generate signatures: one signature per cluster

    Args:
        datasets: List of DiscoveryDatasetSummary objects containing:
            - dataset_id (UUID): Dataset identifier
            - omics_type (str): Type of omics data
            - disease (str): Disease context (optional)
            - matrix (str): Sample type (optional)
            - features (set): Set of feature names in dataset
            - directions (dict): NOT USED in v1 (always None)

        min_support: Minimum number of datasets a feature must appear in to be considered.
            - Lower (e.g., 2): More exploratory, captures rare patterns
            - Higher (e.g., 5): More conservative, only strong patterns
            - Default: 2

        min_overlap: Minimum overlap ratio for clustering features together.
            - Calculated as: len(shared_datasets) / max(len(datasets_for_either))
            - Range: 0.0 to 1.0
            - Lower (e.g., 0.2): Looser clusters, more features per signature
            - Higher (e.g., 0.7): Tighter clusters, fewer features per signature
            - Default: 0.3

    Returns:
        List[DiscoveredSignature]: Candidate signatures, each containing:
            - name (str): Auto-generated name (AUTO_{omics_type}_{hash})
            - modality (str): Omics type
            - components (List[Component]): List of features (direction=None in v1)
            - support (int): Number of datasets where ALL components appear
            - provenance (dict): Source dataset UUIDs

    Example:
        >>> summaries = [
        ...     DiscoveryDatasetSummary(
        ...         dataset_id=uuid1,
        ...         omics_type="lipidomics",
        ...         features={"Cer(d18:1/16:0)", "Cer(d18:1/18:0)", "SM(d18:1/16:0)"}
        ...     ),
        ...     DiscoveryDatasetSummary(
        ...         dataset_id=uuid2,
        ...         omics_type="lipidomics",
        ...         features={"Cer(d18:1/16:0)", "Cer(d18:1/18:0)", "Cer(d18:1/24:0)"}
        ...     )
        ... ]
        >>> candidates = discover_signatures_from_datasets(summaries, min_support=2, min_overlap=0.5)
        >>> for sig in candidates:
        ...     print(f"{sig.name}: {len(sig.components)} components, support={sig.support}")

    Notes:
        - This is v1 (experimental) - outputs require validation
        - Direction analysis not implemented (planned for v2.0)
        - Statistical testing not implemented (planned for v2.0)
        - See docs/SIGNATURE_DISCOVERY.md for validation workflow

    Warning:
        Discovered signatures are CANDIDATES only. Do not use in production without:
        1. Manual inspection for biological plausibility
        2. Literature cross-checking
        3. Statistical validation
        4. Domain expert review
        See docs/SIGNATURE_DISCOVERY.md for complete validation process.
    """
    # STEP 1: Build inverted index (feature → set of datasets containing it)
    # This allows fast lookup of which datasets contain each feature
    feature_to_datasets = defaultdict(set)
    for ds in datasets:
        for feature in ds.features:
            feature_to_datasets[feature].add(ds.dataset_id)

    # STEP 2: Filter features by minimum support
    # Only keep features that appear in at least min_support datasets
    # This removes noise and focuses on recurring patterns
    signature_sets = []
    for feature, dsids in feature_to_datasets.items():
        if len(dsids) >= min_support:
            signature_sets.append((feature, dsids))

    # STEP 3: Cluster features with similar co-occurrence patterns
    # Overlap ratio = shared_datasets / max(datasets_for_either_feature)
    # Features with high overlap (≥ min_overlap) are grouped into clusters
    clusters = []
    seen = set()  # Track features already assigned to clusters

    for i, (feat, ids) in enumerate(signature_sets):
        if feat in seen:
            continue  # Feature already in a cluster

        # Start new cluster with this feature
        cluster = [feat]

        # Find other features that co-occur with this feature
        for j, (other_feat, other_ids) in enumerate(signature_sets):
            if i != j:  # Don't compare feature to itself
                # Calculate overlap ratio
                overlap_ratio = len(ids & other_ids) / max(len(ids), 1)

                if overlap_ratio >= min_overlap:
                    # Features co-occur frequently enough → add to cluster
                    cluster.append(other_feat)
                    seen.add(other_feat)

        seen.add(feat)
        clusters.append(cluster)

    # STEP 4: Generate candidate signatures from clusters
    results = []
    for cluster in clusters:
        # Find datasets where ALL cluster features appear together
        # This is the "support" for the signature
        involved_datasets = (
            set.intersection(*(feature_to_datasets[f] for f in cluster))
            if cluster
            else set()
        )

        # Generate unique signature name using hash of components
        hash_base = " ".join(sorted(cluster)).encode()
        name = f"AUTO_{datasets[0].omics_type}_{hashlib.md5(hash_base).hexdigest()[:8]}"

        # Create components (NOTE: direction=None in v1 - not extracted)
        components = [Component(feature=f) for f in cluster]

        # Create signature candidate
        sig = DiscoveredSignature(
            name=name,
            modality=datasets[0].omics_type,
            components=components,
            support=len(involved_datasets),  # Number of datasets with ALL components
            provenance={"dataset_ids": [str(d) for d in involved_datasets]},
        )
        results.append(sig)

    return results
