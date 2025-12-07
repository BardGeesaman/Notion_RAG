"""
Feature extraction functions.

Extracts features (metabolites, etc.) from mwTab JSON data and text content.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Set

from amprenta_rag.ingestion.features.constants import AMINO_ACIDS, NUCLEOTIDES
from amprenta_rag.ingestion.features.normalization import normalize_metabolite_name


def extract_features_from_mwtab(mwtab_json: Dict[str, Any]) -> List[str]:
    """
    Extract metabolite names from mwTab JSON data.

    Args:
        mwtab_json: Parsed mwTab JSON dictionary

    Returns:
        List of normalized metabolite names
    """
    metabolite_names: Set[str] = set()

    # Check for MS_METABOLITE_DATA section
    metabolite_sections = [
        "MS_METABOLITE_DATA",
        "GC_METABOLITE_DATA",
        "LC_METABOLITE_DATA",
        "METABOLITE_DATA",
    ]

    for section_key in metabolite_sections:
        if section_key not in mwtab_json:
            continue

        section = mwtab_json[section_key]
        if not isinstance(section, dict):
            continue

        # Extract data array
        data_array = section.get("Data", [])
        if not isinstance(data_array, list):
            continue

        for row in data_array:
            if isinstance(row, dict):
                # Look for Metabolite key (case-insensitive)
                metabolite_key = None
                for key in row.keys():
                    if key.lower() in ["metabolite", "metabolite_name", "compound"]:
                        metabolite_key = key
                        break

                if metabolite_key:
                    metabolite_raw = row.get(metabolite_key)
                    if metabolite_raw and isinstance(metabolite_raw, str):
                        normalized = normalize_metabolite_name(metabolite_raw)
                        if normalized:
                            metabolite_names.add(normalized)

            elif isinstance(row, list) and len(row) > 0:
                # First element might be metabolite name
                if isinstance(row[0], str):
                    normalized = normalize_metabolite_name(row[0])
                    if normalized:
                        metabolite_names.add(normalized)

    return sorted(list(metabolite_names))


def extract_features_from_text(text: str) -> List[str]:
    """
    Extract metabolite names from text content using pattern matching.

    This is a lightweight scanner for literature chunks, email/note chunks,
    and experiment descriptions.

    Args:
        text: Text content to scan

    Returns:
        List of normalized metabolite names found in text
    """
    if not text:
        return []

    text_lower = text.lower()
    found_metabolites: Set[str] = set()

    # Scan for common amino acids
    for aa in AMINO_ACIDS:
        # Look for whole word matches (not substring)
        pattern = r"\b" + re.escape(aa) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(aa))

    # Scan for nucleotides
    for nt in NUCLEOTIDES:
        pattern = r"\b" + re.escape(nt) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(nt))

    # Scan for ceramide patterns (basic)
    ceramide_patterns = [
        r"\bcer\s*\([^)]+\)",
        r"\bceramide\s*\([^)]+\)",
        r"\bcer\([^)]+\)",
    ]
    for pattern in ceramide_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            found_metabolites.add(normalize_metabolite_name(match.group(0)))

    # Scan for common small molecules mentioned in signatures
    common_metabolites = [
        "glucose",
        "lactate",
        "pyruvate",
        "citrate",
        "succinate",
        "fumarate",
        "malate",
        "oxaloacetate",
        "alpha-ketoglutarate",
        "acetylate",
        "carnitine",
        "creatine",
        "choline",
    ]
    for metabolite in common_metabolites:
        pattern = r"\b" + re.escape(metabolite) + r"\b"
        if re.search(pattern, text_lower, re.IGNORECASE):
            found_metabolites.add(normalize_metabolite_name(metabolite))

    return sorted(list(found_metabolites))
