"""
Metabolite name normalization for metabolomics ingestion.

This module provides functions for normalizing metabolite names from various
formats into canonical form suitable for knowledge graph linking.
"""

from __future__ import annotations

import re

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def normalize_metabolite_name(raw: str) -> str:
    """
    Normalize metabolite names into a canonical form (best-effort).

    Phase 1: light normalization + synonym cleanup.

    Args:
        raw: Raw metabolite name from file

    Returns:
        Normalized/canonical metabolite name
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Remove adducts: [M+H], [M-H], [M+Na], etc.
    normalized = re.sub(r"\[M[+-][A-Za-z0-9]+\]", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\[M\+H\]\+", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\[M-H\]-", "", normalized, flags=re.IGNORECASE)
    # Remove standalone charge indicators: +, -, +1, -1, etc.
    normalized = re.sub(r"\s*[+-]\d*\s*$", "", normalized)
    normalized = re.sub(r"\s*\+\s*$", "", normalized)
    normalized = re.sub(r"\s*-\s*$", "", normalized)

    # Remove trailing annotations: (pos), (neg), (+), (-)
    normalized = re.sub(r"\s*\(pos\)", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\s*\(neg\)", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\s*\(\+\)", "", normalized)
    normalized = re.sub(r"\s*\(-\)", "", normalized)

    # Replace underscores with spaces
    normalized = normalized.replace("_", " ")

    # Trim again after removals
    normalized = normalized.strip()

    # Basic synonym mapping (optional, minimal for Phase 1)
    synonym_map = {
        "l-glutamic acid": "Glutamate",
        "glutamic acid": "Glutamate",
        "l-glutamine": "Glutamine",
        "glutamine": "Glutamine",
        "l-serine": "Serine",
        "serine": "Serine",
        "l-alanine": "Alanine",
        "alanine": "Alanine",
        "l-aspartic acid": "Aspartate",
        "aspartic acid": "Aspartate",
        "l-aspartate": "Aspartate",
        "aspartate": "Aspartate",
        "l-lysine": "Lysine",
        "lysine": "Lysine",
        "l-arginine": "Arginine",
        "arginine": "Arginine",
        "l-proline": "Proline",
        "proline": "Proline",
        "l-valine": "Valine",
        "valine": "Valine",
        "l-leucine": "Leucine",
        "leucine": "Leucine",
        "l-isoleucine": "Isoleucine",
        "isoleucine": "Isoleucine",
        "l-methionine": "Methionine",
        "methionine": "Methionine",
        "l-phenylalanine": "Phenylalanine",
        "phenylalanine": "Phenylalanine",
        "l-tyrosine": "Tyrosine",
        "tyrosine": "Tyrosine",
        "l-tryptophan": "Tryptophan",
        "tryptophan": "Tryptophan",
        "l-cysteine": "Cysteine",
        "cysteine": "Cysteine",
        "l-threonine": "Threonine",
        "threonine": "Threonine",
        "l-histidine": "Histidine",
        "histidine": "Histidine",
    }

    normalized_lower = normalized.lower()
    if normalized_lower in synonym_map:
        normalized = synonym_map[normalized_lower]

    # Case normalization: Title Case for common metabolites
    # Keep as-is if already looks formatted
    if normalized and not normalized[0].isupper():
        # Convert to title case
        normalized = normalized.title()

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][METABOLOMICS] Normalized metabolite '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][METABOLOMICS] Keeping metabolite name as-is: '%s'",
            original,
        )

    return normalized

