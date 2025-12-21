"""
Lipid species normalization for lipidomics ingestion.

This module provides functions for normalizing raw lipid names from various
vendor formats into canonical species format suitable for knowledge graph linking.
"""

from __future__ import annotations

import re
from typing import Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def normalize_lipid_species(raw_name: str) -> Optional[str]:
    """
    Normalize raw lipid name to canonical species format.

    Extended version of map_raw_lipid_to_canonical_species with support for
    more vendor formats and edge cases.

    Examples:
        "CER 16:0" → "Cer(d18:1/16:0)"
        "Cer 16:0" → "Cer(d18:1/16:0)"
        "Cer d18:1/16:0" → "Cer(d18:1/16:0)"
        "SM 24:1;O2" → "SM(d18:1/24:1)"
        "SM(d18:1_24:1)" → "SM(d18:1/24:1)"
        "hex_cer_24_0" → "HexCer(d18:1/24:0)"
        "CER(d18:1/16:0)+H" → "Cer(d18:1/16:0)"

    Args:
        raw_name: Raw lipid name from file

    Returns:
        Canonical species name or None if cannot be normalized
    """
    if not raw_name or not isinstance(raw_name, str):
        return None

    raw = raw_name.strip()
    if not raw:
        return None

    # If already in canonical format, return as-is (with cleanup)
    if "(" in raw and ")" in raw:
        # Clean up: remove adducts, fix separators
        canonical = re.sub(r"[\+\-].*$", "", raw)  # Remove adducts (+H, -H2O, etc.)
        canonical = canonical.replace("_", "/")  # Fix underscore separators
        canonical = canonical.replace(" ", "")  # Remove spaces
        # Ensure proper format: Class(d18:1/chain) or Class(d18:1_chain) -> Class(d18:1/chain)
        if re.match(r"^[A-Za-z]+\(d\d+:\d+[/_]\d+:\d+\)$", canonical):
            # Normalize separator to /
            canonical = canonical.replace("_", "/")
            logger.debug(
                "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
                raw_name,
                canonical,
            )
            return canonical

    raw_lower = raw.lower().strip()

    # Remove adducts and charge info
    raw_lower = re.sub(r"[\+\-].*$", "", raw_lower)
    raw_lower = re.sub(r";.*$", "", raw_lower)  # Remove modifications like ;O2

    # Class mapping
    class_map = {
        "cer": "Cer",
        "ceramide": "Cer",
        "sm": "SM",
        "sphingomyelin": "SM",
        "hexcer": "HexCer",
        "hexosylceramide": "HexCer",
        "laccer": "LacCer",
        "lactosylceramide": "LacCer",
        "glccer": "GlcCer",
        "glucosylceramide": "GlcCer",
    }

    result: Optional[str] = None

    # Try to match various formats
    # Format 1: CER 16:0, SM 24:1, etc.
    simple_match = re.match(
        r"([a-z]+)\s*([cC]?)(\d+):(\d+)", raw_lower
    )
    if simple_match:
        class_name = simple_match.group(1)
        chain = f"{simple_match.group(3)}:{simple_match.group(4)}"
        canonical_class = class_map.get(class_name, class_name.capitalize())
        result = f"{canonical_class}(d18:1/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Format 2: Cer d18:1/16:0, SM d18:1/24:1, etc.
    with_backbone = re.match(
        r"([a-z]+)\s*d(\d+):(\d+)[_/](\d+):(\d+)", raw_lower
    )
    if with_backbone:
        class_name = with_backbone.group(1)
        backbone = f"d{with_backbone.group(2)}:{with_backbone.group(3)}"
        chain = f"{with_backbone.group(4)}:{with_backbone.group(5)}"
        canonical_class = class_map.get(class_name, class_name.capitalize())
        result = f"{canonical_class}({backbone}/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Format 3: hex_cer_24_0, sm_18_1, etc.
    underscore_format = re.match(
        r"([a-z_]+)_(\d+)_(\d+)", raw_lower
    )
    if underscore_format:
        class_part = underscore_format.group(1)
        chain = f"{underscore_format.group(2)}:{underscore_format.group(3)}"
        # Extract class name - check for hex, lac, glc prefixes first
        if "hex" in class_part:
            result = f"HexCer(d18:1/{chain})"
        elif "lac" in class_part:
            result = f"LacCer(d18:1/{chain})"
        elif "glc" in class_part:
            result = f"GlcCer(d18:1/{chain})"
        else:
            # Check other class mappings
            for key, canonical_class in class_map.items():
                if key in class_part:
                    result = f"{canonical_class}(d18:1/{chain})"
                    break
            else:
                # Default to Cer if no match
                result = f"Cer(d18:1/{chain})"
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    # Try the existing map_raw_lipid_to_canonical_species as fallback
    from amprenta_rag.ingestion.signature_matching import (
        map_raw_lipid_to_canonical_species)
    result = map_raw_lipid_to_canonical_species(raw_name)
    if result is not None:
        logger.info(
            "[INGEST][LIPIDOMICS] Normalized raw lipid '%s' -> '%s'",
            raw_name,
            result,
        )
        return result

    return raw_name

