"""
File parsing utilities for proteomics ingestion.

This module provides functions for parsing CSV/TSV files and extracting
protein/gene identifiers with automatic column detection.
"""

from __future__ import annotations

from typing import Set, Tuple

import pandas as pd

from amprenta_rag.ingestion.proteomics.normalization import normalize_protein_identifier
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_protein_set_from_file(file_path: str) -> Tuple[Set[str], int]:
    """
    Load CSV/TSV and return a set of normalized protein/gene identifiers.

    Args:
        file_path: Path to the proteomics file

    Returns:
        Tuple of (set of normalized proteins, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][PROTEOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][PROTEOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][PROTEOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect protein/gene identity column
    protein_column = None
    candidate_names = [
        "Protein",
        "protein",
        "Protein ID",
        "ProteinID",
        "protein_id",
        "Gene",
        "gene",
        "Gene Symbol",
        "GeneSymbol",
        "gene_name",
        "protein_name",
        "Protein Name",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            protein_column = col
            break

    if not protein_column:
        # Try case-insensitive match
        for col in df.columns:
            col_lower = col.lower()
            if any(
                cand.lower() in col_lower or col_lower in cand.lower()
                for cand in candidate_names
            ):
                protein_column = col
                break

    if not protein_column:
        logger.error(
            "[INGEST][PROTEOMICS] Could not find protein/gene identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find protein/gene identity column in {file_path}"
        )

    logger.info(
        "[INGEST][PROTEOMICS] Using column '%s' for protein identity in file %s",
        protein_column,
        file_path,
    )

    # Extract and normalize proteins
    protein_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(protein_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_protein_identifier(raw_name_str)
        if normalized:
            protein_set.add(normalized)

    logger.info(
        "[INGEST][PROTEOMICS] Extracted %d unique proteins from %d rows in file %s",
        len(protein_set),
        total_rows,
        file_path,
    )

    return protein_set, total_rows

