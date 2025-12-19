"""
Gene identifier normalization for transcriptomics ingestion.

This module provides functions for normalizing gene identifiers from various
formats into canonical gene symbol format suitable for knowledge graph linking.
"""

from __future__ import annotations

import re

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def normalize_gene_identifier(raw: str) -> str:
    """
    Normalize gene identifiers into a canonical gene symbol-like format.

    Phase 1: best-effort cleaning to HGNC-style symbols where possible.

    Handles:
    - Species suffixes: TP53_HUMAN → TP53
    - Case normalization: tp53 → TP53
    - Annotations: TP53 (Human) → TP53
    - Ensembl IDs: Keep as-is but clean whitespace

    Args:
        raw: Raw gene identifier from file

    Returns:
        Normalized/canonical gene identifier
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Check if it's an Ensembl ID (ENSG0000... or similar)
    ensembl_pattern = r"^ENS[A-Z]*G\d+"
    is_ensembl = bool(re.match(ensembl_pattern, normalized, re.IGNORECASE))

    if is_ensembl:
        # Keep Ensembl IDs as-is but ensure uppercase and clean whitespace
        normalized = normalized.upper().strip()
        if normalized != original:
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Normalized Ensembl ID '%s' -> '%s'",
                original,
                normalized,
            )
        else:
            logger.debug(
                "[INGEST][TRANSCRIPTOMICS] Keeping Ensembl ID as-is: '%s'",
                original,
            )
        return normalized

    # Remove species suffixes: TP53_HUMAN, Tp53_mouse → TP53
    # Common species suffixes: HUMAN, MOUSE, RAT, BOVIN, PIG, CHICK
    if "_" in normalized:
        parts = normalized.split("_")
        if len(parts) > 1:
            last_part = parts[-1].upper()
            if last_part in ["HUMAN", "MOUSE", "RAT", "BOVIN", "PIG", "CHICK", "HSA", "MMU"]:
                normalized = "_".join(parts[:-1])

    # Remove bracketed annotations: TP53 (Human) → TP53
    normalized = re.sub(r"\s*\([^)]+\)", "", normalized)
    normalized = re.sub(r"\s*\[[^\]]+\]", "", normalized)

    # Remove trailing version-like suffixes if clearly non-gene: TP53.1 → TP53
    # But be careful not to remove legitimate gene symbols that end in numbers
    # Only remove if it looks like a version suffix (e.g., .1, .2, .v1)
    if re.match(r"^[A-Z]+\d*\.\d+$", normalized, re.IGNORECASE):
        normalized = re.sub(r"\.\d+$", "", normalized)

    # Convert to uppercase for gene symbols (tp53 → TP53)
    # Gene symbols are typically uppercase letters and numbers
    if re.match(r"^[A-Za-z]+\d*[A-Za-z]*$", normalized):
        normalized = normalized.upper()

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Normalized '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][TRANSCRIPTOMICS] Keeping gene identifier as-is: '%s'",
            original,
        )

    return normalized

