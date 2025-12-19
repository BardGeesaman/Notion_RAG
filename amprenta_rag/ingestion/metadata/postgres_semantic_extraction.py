"""
Postgres-compatible semantic metadata extraction.

Extracts semantic metadata (diseases, targets, lipid signatures, etc.) from text content
without requiring Notion pages. This is a Postgres-only alternative to the Notion-based
semantic metadata extraction.

This module provides basic pattern-matching extraction and can be enhanced with
LLM-based extraction for better accuracy. Set ENABLE_LLM_SEMANTIC_EXTRACTION=true to enable.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Try to import LLM extraction (optional enhancement)
try:
    pass

    LLM_AVAILABLE = True
except ImportError:
    LLM_AVAILABLE = False


# Common disease terms for pattern matching
DISEASE_PATTERNS = {
    "Fragile X syndrome": [r"\bfxs\b", r"\bfragile x\b", r"\bfragile-x\b", r"\bfxs\b"],
    "ALS": [r"\bals\b", r"\bamyotrophic\b", r"\bamyotrophic lateral sclerosis\b"],
    "Alzheimer's disease": [r"\balzheimer\b", r"\bad\b", r"\balzheimer's\b"],
    "Parkinson's disease": [r"\bparkinson\b", r"\bpd\b", r"\bparkinson's\b"],
    "Huntington's disease": [r"\bhuntington\b", r"\bhd\b", r"\bhuntington's\b"],
}

# Common target patterns
TARGET_PATTERNS = [
    r"\b(\w+)\s*(?:protein|receptor|enzyme|kinase|inhibitor)\b",
    r"\b(?:target|drug target|therapeutic target):\s*(\w+)",
]

# Common signature patterns
SIGNATURE_PATTERNS = [
    r"\b(lipid signature|metabolic signature|biomarker signature|signature)\s*:?\s*(\w+)",
    r"\bsignature\s+(?:id|short id|name):\s*(\w+)",
]


def extract_diseases_from_text(text: str, case_sensitive: bool = False) -> List[str]:
    """
    Extract disease terms from text using pattern matching.

    Args:
        text: Text content to search
        case_sensitive: Whether to match case-sensitively

    Returns:
        List of disease names found
    """
    if not text:
        return []

    search_text = text if case_sensitive else text.lower()
    diseases_found = []

    for disease_name, patterns in DISEASE_PATTERNS.items():
        for pattern in patterns:
            if case_sensitive:
                if re.search(pattern, text, re.IGNORECASE):
                    if disease_name not in diseases_found:
                        diseases_found.append(disease_name)
                    break
            else:
                if re.search(pattern, search_text):
                    if disease_name not in diseases_found:
                        diseases_found.append(disease_name)
                    break

    return diseases_found


def extract_targets_from_text(text: str) -> List[str]:
    """
    Extract target molecules/proteins from text using pattern matching.

    Args:
        text: Text content to search

    Returns:
        List of target names found
    """
    if not text:
        return []

    targets_found = []

    for pattern in TARGET_PATTERNS:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            target = match.group(1) if match.groups() else match.group(0)
            if target and target not in targets_found:
                # Clean up the target name
                target = target.strip().title()
                if len(target) > 2:  # Filter out very short matches
                    targets_found.append(target)

    return targets_found


def extract_signatures_from_text(text: str) -> List[str]:
    """
    Extract signature identifiers from text using pattern matching.

    Args:
        text: Text content to search

    Returns:
        List of signature identifiers found
    """
    if not text:
        return []

    signatures_found = []

    for pattern in SIGNATURE_PATTERNS:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            sig = match.group(2) if len(match.groups()) >= 2 else match.group(1)
            if sig and sig not in signatures_found:
                signatures_found.append(sig.strip())

    return signatures_found


def extract_semantic_metadata_from_text(
    text: str,
    extract_diseases: bool = True,
    extract_targets: bool = True,
    extract_signatures: bool = True,
) -> Dict[str, Any]:
    """
    Extract semantic metadata from text content (Postgres-compatible).

    This function extracts structured metadata from unstructured text without
    requiring Notion pages. It uses pattern matching and can be enhanced with
    LLM-based extraction for better accuracy.

    Args:
        text: Text content to extract metadata from
        extract_diseases: Whether to extract disease terms
        extract_targets: Whether to extract target molecules
        extract_signatures: Whether to extract signature identifiers

    Returns:
        Dictionary with semantic metadata:
        - diseases: List[str]
        - targets: List[str]
        - lipid_signatures: List[str] (if signatures found)
    """
    metadata: Dict[str, Any] = {
        "diseases": [],
        "targets": [],
        "lipid_signatures": [],
    }

    if not text:
        return metadata

    # Extract diseases
    if extract_diseases:
        diseases = extract_diseases_from_text(text)
        if diseases:
            metadata["diseases"] = diseases
            logger.debug(
                "[SEMANTIC] Extracted %d disease(s) from text: %s",
                len(diseases),
                diseases,
            )

    # Extract targets
    if extract_targets:
        targets = extract_targets_from_text(text)
        if targets:
            metadata["targets"] = targets
            logger.debug(
                "[SEMANTIC] Extracted %d target(s) from text: %s",
                len(targets),
                targets,
            )

    # Extract signatures
    if extract_signatures:
        signatures = extract_signatures_from_text(text)
        if signatures:
            metadata["lipid_signatures"] = signatures
            logger.debug(
                "[SEMANTIC] Extracted %d signature(s) from text: %s",
                len(signatures),
                signatures,
            )

    return metadata


def merge_semantic_metadata(
    existing: Dict[str, Any],
    extracted: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Merge extracted semantic metadata with existing metadata.

    Args:
        existing: Existing metadata dictionary (from Postgres fields)
        extracted: Newly extracted metadata from text

    Returns:
        Merged metadata dictionary
    """
    merged = existing.copy()

    # Merge diseases
    existing_diseases = set(merged.get("diseases", []))
    extracted_diseases = set(extracted.get("diseases", []))
    merged["diseases"] = sorted(list(existing_diseases | extracted_diseases))

    # Merge targets
    existing_targets = set(merged.get("targets", []))
    extracted_targets = set(extracted.get("targets", []))
    merged["targets"] = sorted(list(existing_targets | extracted_targets))

    # Merge signatures
    existing_sigs = set(merged.get("lipid_signatures", []))
    extracted_sigs = set(extracted.get("lipid_signatures", []))
    merged["lipid_signatures"] = sorted(list(existing_sigs | extracted_sigs))

    return merged


def enhance_metadata_with_semantic_extraction(
    metadata: Dict[str, Any],
    text_content: str,
    source_type: str = "generic",
    use_llm: bool = False,
) -> Dict[str, Any]:
    """
    Enhance existing metadata by extracting semantic information from text.

    This is a convenience function that extracts semantic metadata from text
    and merges it with existing metadata. Can optionally use LLM for more accurate extraction.

    Args:
        metadata: Existing metadata dictionary
        text_content: Text content to extract from
        source_type: Type of source (dataset, experiment, email, literature) - for LLM logging
        use_llm: Whether to use LLM-based extraction (if available and enabled)

    Returns:
        Enhanced metadata dictionary
    """
    from amprenta_rag.config import get_config

    cfg = get_config()

    # Try LLM extraction if requested and available
    if use_llm and LLM_AVAILABLE and cfg.pipeline.enable_llm_semantic_extraction:
        try:
            from amprenta_rag.ingestion.metadata.llm_semantic_extraction import (
                enhance_metadata_with_llm,
            )
            return enhance_metadata_with_llm(metadata, text_content, source_type)
        except Exception as e:
            logger.debug(
                "[SEMANTIC] LLM extraction failed, falling back to pattern matching: %r",
                e,
            )

    # Fall back to pattern matching
    extracted = extract_semantic_metadata_from_text(text_content)
    return merge_semantic_metadata(metadata, extracted)

