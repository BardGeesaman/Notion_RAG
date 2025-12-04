"""
Protein identifier normalization for proteomics ingestion.

This module provides functions for normalizing protein/gene identifiers from various
formats into canonical format suitable for knowledge graph linking.
"""

from __future__ import annotations

import re

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def normalize_protein_identifier(raw: str) -> str:
    """
    Normalize protein or gene identifiers into a canonical symbol-like format.

    Phase 1: basic mapping + cleaning.

    Handles:
    - FASTA prefixes: sp|P12345|TP53_HUMAN → TP53
    - Isoform suffixes: Q9Y6K9-2 → Q9Y6K9
    - Gene symbols: actb → ACTB
    - Annotations: TP53 (Human) → TP53

    Args:
        raw: Raw protein/gene identifier from file

    Returns:
        Normalized/canonical identifier
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Remove FASTA prefixes: sp|P12345|TP53_HUMAN → TP53_HUMAN
    # Pattern: sp|UNIPROT_ID|GENE_SPECIES or tr|...|...
    fasta_match = re.match(r"^(sp|tr|ref)\|[^|]+\|([^|]+)", normalized)
    if fasta_match:
        normalized = fasta_match.group(2)
        logger.debug(
            "[INGEST][PROTEOMICS] Extracted from FASTA format: '%s' -> '%s'",
            original,
            normalized,
        )

    # Extract gene symbol from format like TP53_HUMAN → TP53
    # Or P04637_TP53_HUMAN → TP53
    if "_" in normalized:
        parts = normalized.split("_")
        # If last part looks like species (HUMAN, MOUSE, etc.), remove it
        if len(parts) > 1 and parts[-1].upper() in [
            "HUMAN",
            "MOUSE",
            "RAT",
            "BOVIN",
            "PIG",
            "CHICK",
        ]:
            normalized = "_".join(parts[:-1])
        # If we have something like P04637_TP53, try to extract TP53
        if len(parts) > 1:
            # Check if first part looks like UniProt ID (alphanumeric, 6-10 chars)
            if re.match(r"^[A-Z0-9]{6,10}$", parts[0].upper()):
                # Use the second part as gene symbol
                if len(parts) > 1:
                    normalized = parts[1]

    # Remove isoform suffixes: Q9Y6K9-2 → Q9Y6K9
    # But keep if it's part of a gene symbol pattern
    if re.match(r"^[A-Z0-9]+-\d+$", normalized.upper()):
        # Looks like UniProt ID with isoform, remove suffix
        normalized = re.sub(r"-\d+$", "", normalized)
    elif "-" in normalized and not re.match(r"^[A-Z]+\d+[A-Z]*$", normalized.upper()):
        # Has hyphen but not a simple gene symbol, might be isoform
        # Only remove if it's at the end and looks like a number
        normalized = re.sub(r"-\d+$", "", normalized)

    # Remove bracketed annotations: TP53 (Human) → TP53
    normalized = re.sub(r"\s*\([^)]+\)", "", normalized)
    normalized = re.sub(r"\s*\[[^\]]+\]", "", normalized)

    # Convert to uppercase for gene symbols (if it looks like a gene symbol)
    # Gene symbols are typically 1-10 uppercase letters/numbers
    if re.match(r"^[A-Za-z0-9]{1,10}$", normalized):
        normalized = normalized.upper()
    # If it's a UniProt ID (6-10 alphanumeric), keep uppercase
    elif re.match(r"^[A-Z0-9]{6,10}$", normalized.upper()):
        normalized = normalized.upper()

    # Replace underscores with hyphens only if it looks gene-like
    # (Not for UniProt IDs which may have underscores)
    if re.match(r"^[A-Z]+\d*[A-Z]*$", normalized.upper()):
        normalized = normalized.replace("_", "-")

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][PROTEOMICS] Normalized '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][PROTEOMICS] Keeping identifier as-is: '%s'",
            original,
        )

    return normalized

