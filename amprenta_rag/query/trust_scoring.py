"""Hierarchical RAG trust scoring for source reliability."""
from __future__ import annotations

from typing import Dict, List, Any

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

TRUST_LEVELS = {
    "internal_validated": 1.0,  # Curated internal data
    "internal": 0.8,            # Internal experiments/results
    "external_curated": 0.6,    # Trusted external (GEO, MetaboLights)
    "external": 0.4,            # Other external sources
    "literature": 0.5,          # Published papers
    "general": 0.2,             # General knowledge
}


def get_source_trust(source_type: str, metadata: Dict[str, Any] = None) -> float:
    """
    Get trust score for a source based on type and metadata.

    Args:
        source_type: Type of source (e.g., "Experiment", "Dataset", "Literature")
        metadata: Optional metadata dict with validation flags

    Returns:
        Trust score between 0.0 and 1.0
    """
    if metadata is None:
        metadata = {}

    # Normalize source_type to lowercase for matching
    source_lower = (source_type or "").lower()

    # Check for validation flags in metadata that boost trust
    is_validated = metadata.get("validated", False) or metadata.get("qc_passed", False)
    is_curated = metadata.get("curated", False) or metadata.get("reviewed", False)
    is_internal = metadata.get("internal", False) or source_lower in ["experiment", "dataset", "signature"]

    # Determine base trust level
    if is_internal and is_validated:
        base_trust = TRUST_LEVELS["internal_validated"]
    elif is_internal:
        base_trust = TRUST_LEVELS["internal"]
    elif is_curated and source_lower in ["geo", "metabolights", "arrayexpress", "pride"]:
        base_trust = TRUST_LEVELS["external_curated"]
    elif source_lower in ["literature", "paper", "publication", "zotero"]:
        base_trust = TRUST_LEVELS["literature"]
    elif source_lower in ["email", "note", "general"]:
        base_trust = TRUST_LEVELS["general"]
    else:
        base_trust = TRUST_LEVELS["external"]

    # Apply boosts from metadata
    if is_validated:
        base_trust = min(1.0, base_trust + 0.1)
    if is_curated:
        base_trust = min(1.0, base_trust + 0.05)

    # Check for quality metrics that affect trust
    quality_score = metadata.get("quality_score")
    if quality_score is not None:
        # Normalize quality_score (assuming 0-1 range) and blend with base trust
        base_trust = (base_trust * 0.7) + (quality_score * 0.3)

    return min(1.0, max(0.0, base_trust))


def weight_results_by_trust(matches: List[Dict[str, Any]], alpha: float = 0.3) -> List[Dict[str, Any]]:
    """
    Weight RAG results by trust scores and reorder.

    Args:
        matches: List of match dicts with 'score' and 'metadata' keys
        alpha: Weight for trust score (0.0-1.0). Higher alpha = more weight on trust.
               Default 0.3 means 70% original score, 30% trust score.

    Returns:
        List of matches with 'trust_score' added, reordered by final_score
    """
    if not matches:
        return []

    # Calculate trust scores and combine with original scores
    weighted_matches = []

    for match in matches:
        # Extract score and metadata
        original_score = match.get("score", 0.0)
        metadata = match.get("metadata", {})
        source_type = metadata.get("source_type") or metadata.get("source") or "unknown"

        # Get trust score
        trust_score = get_source_trust(source_type, metadata)

        # Combine scores: final = (1-alpha)*original + alpha*trust
        final_score = (1.0 - alpha) * original_score + alpha * trust_score

        # Create new match dict with trust info
        weighted_match = match.copy()
        weighted_match["trust_score"] = trust_score
        weighted_match["final_score"] = final_score
        weighted_match["original_score"] = original_score

        weighted_matches.append(weighted_match)

    # Sort by final_score descending
    weighted_matches.sort(key=lambda m: m.get("final_score", 0.0), reverse=True)

    return weighted_matches


def get_trust_summary(matches: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Get summary statistics about trust levels in matches.

    Args:
        matches: List of match dicts (should have trust_score if weighted)

    Returns:
        Dictionary with:
            - "levels": dict mapping trust level names to counts
            - "average_trust": average trust score
            - "high_trust_count": number of matches with trust >= 0.7
    """
    if not matches:
        return {
            "levels": {},
            "average_trust": 0.0,
            "high_trust_count": 0,
        }

    # Calculate trust scores if not present
    trust_scores = []
    for match in matches:
        if "trust_score" in match:
            trust_scores.append(match["trust_score"])
        else:
            # Calculate trust score if not present
            metadata = match.get("metadata", {})
            source_type = metadata.get("source_type") or metadata.get("source") or "unknown"
            trust_score = get_source_trust(source_type, metadata)
            trust_scores.append(trust_score)

    # Categorize by trust level
    level_counts = {
        "internal_validated": 0,
        "internal": 0,
        "external_curated": 0,
        "external": 0,
        "literature": 0,
        "general": 0,
    }

    for trust_score in trust_scores:
        if trust_score >= 0.9:
            level_counts["internal_validated"] += 1
        elif trust_score >= 0.75:
            level_counts["internal"] += 1
        elif trust_score >= 0.55:
            level_counts["external_curated"] += 1
        elif trust_score >= 0.45:
            level_counts["literature"] += 1
        elif trust_score >= 0.3:
            level_counts["external"] += 1
        else:
            level_counts["general"] += 1

    # Calculate statistics
    average_trust = sum(trust_scores) / len(trust_scores) if trust_scores else 0.0
    high_trust_count = sum(1 for ts in trust_scores if ts >= 0.7)

    return {
        "levels": level_counts,
        "average_trust": average_trust,
        "high_trust_count": high_trust_count,
    }

