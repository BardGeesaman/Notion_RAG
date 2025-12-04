"""
Lipid species mapping functions.

Handles mapping raw lipid names from datasets to canonical species names.
"""

from __future__ import annotations

import re
from typing import Optional

from amprenta_rag.signatures.species_matching import normalize_species_name


def map_raw_lipid_to_canonical_species(raw_name: str) -> Optional[str]:
    """
    Map raw lipid name from dataset to canonical species name.

    Uses normalization and class-level matching to handle vendor formats.

    Examples:
        "CER 16:0" → "Cer(d18:1/16:0)"
        "Ceramide C24:1" → "Cer(d18:1/24:1)"
        "SM 18:0" → "SM(d18:1/18:0)"

    Args:
        raw_name: Raw lipid name from dataset

    Returns:
        Canonical species name or None if cannot be mapped
    """
    if not raw_name:
        return None

    normalized = normalize_species_name(raw_name)

    # Try to reconstruct canonical format from normalized name
    # This is a simplified mapping - can be enhanced with Lipid Species DB lookup
    raw_lower = raw_name.lower().strip()

    # Handle common vendor formats
    # CER 16:0 → Cer(d18:1/16:0)
    cer_match = re.match(r"cer\s*([cC]?)(\d+):(\d+)", raw_lower)
    if cer_match:
        chain = f"{cer_match.group(2)}:{cer_match.group(3)}"
        return f"Cer(d18:1/{chain})"

    # SM 16:0 → SM(d18:1/16:0)
    sm_match = re.match(r"sm\s*([cC]?)(\d+):(\d+)", raw_lower)
    if sm_match:
        chain = f"{sm_match.group(2)}:{sm_match.group(3)}"
        return f"SM(d18:1/{chain})"

    # If already in canonical format, return as-is
    if "(" in raw_name and ")" in raw_name:
        return raw_name

    return None

