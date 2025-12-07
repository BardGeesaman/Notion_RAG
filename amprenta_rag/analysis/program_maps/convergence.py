"""
Convergence indicator computation for program-signature mapping.

Functions for analyzing cross-omics convergence metrics.
"""

from typing import Dict, List

from amprenta_rag.analysis.program_maps.models import ProgramSignatureScore
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def compute_convergence_indicators(
    signature_scores: List[ProgramSignatureScore],
) -> Dict[str, float]:
    """
    Compute cross-omics convergence indicators.

    Analyzes signature scores to determine:
    - Fraction of signatures with multi-omics support
    - Average number of omics types per signature
    - Count of multi-omics signatures

    Args:
        signature_scores: List of ProgramSignatureScore objects

    Returns:
        Dictionary of convergence metrics:
        - convergence_fraction: Fraction of signatures with >1 omics type
        - avg_omics_per_signature: Average omics types per signature
        - multi_omics_signature_count: Number of multi-omics signatures

    Example:
        >>> indicators = compute_convergence_indicators([score1, score2])
        >>> "convergence_fraction" in indicators
        True
    """
    if not signature_scores:
        return {}

    # Count signatures with multi-omics support
    multi_omics_signatures = 0
    total_omics_types = 0

    for score in signature_scores:
        omics_count = len([s for s in score.score_by_omics.values() if s > 0])
        if omics_count > 1:
            multi_omics_signatures += 1
        total_omics_types += omics_count

    convergence_fraction = multi_omics_signatures / len(signature_scores) if signature_scores else 0.0
    avg_omics_per_signature = total_omics_types / len(signature_scores) if signature_scores else 0.0

    return {
        "convergence_fraction": convergence_fraction,
        "avg_omics_per_signature": avg_omics_per_signature,
        "multi_omics_signature_count": multi_omics_signatures,
    }
