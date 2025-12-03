"""
Species matching utilities for lipid signatures.

Handles matching lipid species between datasets and signatures with:
- Exact matching
- Name normalization
- Class-level fallback
- Optional RefMet mapping
"""

import re
from collections import defaultdict
from typing import Dict, Optional, Set


def normalize_species_name(name: str) -> str:
    """
    Normalize a lipid species name for matching.

    Examples:
        Cer(d18:1/24:1) → ceramide_d18_1_24_1
        Cer C24:1 → ceramide_c24_1
        SM(d18:1/16:0) → sphingomyelin_d18_1_16_0

    Args:
        name: Original species name

    Returns:
        Normalized name string
    """
    if not name:
        return ""

    # Convert to lowercase
    normalized = name.lower().strip()

    # Handle common patterns
    # Cer(d18:1/24:1) → ceramide_d18_1_24_1
    normalized = re.sub(r"cer\(([^)]+)\)", r"ceramide_\1", normalized)
    normalized = normalized.replace("(", "_").replace(")", "")
    normalized = normalized.replace("/", "_").replace(":", "_")
    normalized = normalized.replace("-", "_").replace(" ", "_")

    # Handle "Cer C24:1" → "ceramide_c24_1"
    normalized = re.sub(r"\bcer\s+c(\d+):(\d+)\b", r"ceramide_c\1_\2", normalized)
    normalized = re.sub(r"\bcer\s+(\d+):(\d+)\b", r"ceramide_\1_\2", normalized)

    # Handle SM patterns
    normalized = re.sub(r"\bsm\(([^)]+)\)", r"sphingomyelin_\1", normalized)
    normalized = re.sub(r"\bsm\s+", "sphingomyelin_", normalized)

    # Handle HexCer, LacCer patterns
    normalized = re.sub(r"\bhexcer\(([^)]+)\)", r"hexosylceramide_\1", normalized)
    normalized = re.sub(r"\blacc?er\(([^)]+)\)", r"lactosylceramide_\1", normalized)

    # Clean up multiple underscores
    normalized = re.sub(r"_+", "_", normalized)
    normalized = normalized.strip("_")

    return normalized


def classify_lipid_class(name: str) -> Optional[str]:
    """
    Classify a lipid species into a lipid class.

    Args:
        name: Lipid species name

    Returns:
        Lipid class name or None if not recognized
    """
    if not name:
        return None

    lower = name.lower()

    # Order matters: more specific first
    if "ceramide-1-phosphate" in lower or "c1p" in lower:
        return "Ceramide-1-phosphate"
    if "dihydroceramide" in lower or "dhcer" in lower:
        return "Dihydroceramide"
    if "hexcer" in lower or "glccer" in lower or "galcer" in lower:
        return "Hexosylceramide"
    if "lacc" in lower or "laccer" in lower or "lactosylceramide" in lower:
        return "Lactosylceramide"
    if "sphingomyelin" in lower or lower.startswith("sm(") or lower.startswith("sm "):
        return "Sphingomyelin"
    if "sphingosine-1-phosphate" in lower or "s1p" in lower:
        return "Sphingosine-1-phosphate"
    if "sphingosine" in lower and "phosphate" not in lower:
        return "Sphingosine"
    if "deoxysphingolipid" in lower or "deoxycer" in lower:
        return "Deoxysphingolipid"
    if any(
        g in lower for g in ["gm1", "gm2", "gm3", "gd1", "gd2", "gd3", "gt1", "gt2"]
    ):
        return "Ganglioside"
    if "ceramide" in lower or "cer(d" in lower or lower.startswith("cer "):
        return "Ceramide"

    return None


def extract_fatty_acid_chain(name: str) -> Optional[str]:
    """
    Extract fatty acid chain information from a lipid name.

    Examples:
        Cer(d18:1/24:1) → "24:1"
        Cer C24:1 → "24:1"
        SM(d18:1/16:0) → "16:0"

    Args:
        name: Lipid species name

    Returns:
        Fatty acid chain string or None
    """
    if not name:
        return None

    # Look for patterns like /24:1, C24:1, 24:1
    patterns = [
        r"/(\d+:\d+)",  # /24:1
        r"[cC](\d+:\d+)",  # C24:1
        r"(\d+:\d+)",  # 24:1 (fallback)
    ]

    for pattern in patterns:
        match = re.search(pattern, name)
        if match:
            return match.group(1)

    return None


def match_species(
    dataset_species: Set[str],
    signature_species: Set[str],
    refmet_map: Optional[Dict[str, str]] = None,
) -> Dict[str, Optional[str]]:
    """
    Match signature species to dataset species.

    Returns a mapping: signature_species → matched_dataset_species or None

    Matching strategy:
    1. Exact normalized match
    2. RefMet mapping (if available)
    3. Class-level + fatty acid chain match
    4. Class-level match only

    Args:
        dataset_species: Set of species names from the dataset
        signature_species: Set of species names from the signature
        refmet_map: Optional dict mapping RefMet IDs/names to normalized names

    Returns:
        Dictionary mapping signature_species → matched_dataset_species or None
    """
    # Normalize all names
    dataset_normalized = {normalize_species_name(s): s for s in dataset_species}
    signature_normalized = {normalize_species_name(s): s for s in signature_species}

    # Build class and chain mappings for fallback
    dataset_by_class: Dict[str, Dict[str, str]] = defaultdict(dict)
    dataset_by_class_chain: Dict[str, Dict[str, str]] = defaultdict(dict)

    for norm_name, orig_name in dataset_normalized.items():
        lipid_class = classify_lipid_class(orig_name)
        if lipid_class:
            dataset_by_class[lipid_class][norm_name] = orig_name
            chain = extract_fatty_acid_chain(orig_name)
            if chain:
                dataset_by_class_chain[lipid_class][chain] = orig_name

    # Match results
    matches: Dict[str, Optional[str]] = {}

    for sig_norm, sig_orig in signature_normalized.items():
        matched = None

        # 1. Exact normalized match
        if sig_norm in dataset_normalized:
            matched = dataset_normalized[sig_norm]

        # 2. RefMet mapping (if available)
        if not matched and refmet_map:
            for refmet_key, mapped_name in refmet_map.items():
                if normalize_species_name(refmet_key) == sig_norm:
                    mapped_norm = normalize_species_name(mapped_name)
                    if mapped_norm in dataset_normalized:
                        matched = dataset_normalized[mapped_norm]
                        break

        # 3. Class-level + fatty acid chain match
        if not matched:
            sig_class = classify_lipid_class(sig_orig)
            sig_chain = extract_fatty_acid_chain(sig_orig)

            if sig_class and sig_chain and sig_class in dataset_by_class_chain:
                if sig_chain in dataset_by_class_chain[sig_class]:
                    matched = dataset_by_class_chain[sig_class][sig_chain]

        # 4. Class-level match only (last resort)
        if not matched and sig_class and sig_class in dataset_by_class:
            # Use first available species in the class
            matched = next(iter(dataset_by_class[sig_class].values()))

        matches[sig_orig] = matched

    return matches
