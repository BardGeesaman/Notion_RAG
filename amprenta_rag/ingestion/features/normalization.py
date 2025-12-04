"""
Feature name normalization functions.

Handles normalization of metabolite names and other features from
various formats (mwTab, CSV, vendor formats) into canonical forms.
"""

from __future__ import annotations

import re

from amprenta_rag.ingestion.features.constants import METABOLITE_SYNONYMS


def normalize_metabolite_name(raw: str) -> str:
    """
    Normalize metabolite names from mwTab and CSV formats.

    Args:
        raw: Raw metabolite name (e.g., "HMDB:12345 Glutamine", "L-glutamate")

    Returns:
        Canonical normalized name (e.g., "glutamine", "glutamate")
    """
    if not raw:
        return ""

    # Lowercase and strip
    normalized = raw.lower().strip()

    # Remove prefixes like "HMDB:", "KEGG:", "CHEBI:", etc.
    normalized = re.sub(
        r"^(hmdb|kegg|chebi|pubchem|cas)[:\s]+", "", normalized, flags=re.IGNORECASE
    )

    # Remove common prefixes/suffixes
    normalized = re.sub(r"^\d+\s*[-:]?\s*", "", normalized)  # Leading numbers
    normalized = re.sub(r"\s*\(.*?\)\s*$", "", normalized)  # Trailing parentheses

    # Strip again after cleanup
    normalized = normalized.strip()

    # Apply synonym mapping
    normalized_lower = normalized.lower()
    if normalized_lower in METABOLITE_SYNONYMS:
        normalized = METABOLITE_SYNONYMS[normalized_lower]

    # Capitalize first letter for canonical form
    if normalized:
        normalized = (
            normalized[0].upper() + normalized[1:]
            if len(normalized) > 1
            else normalized.upper()
        )

    return normalized

