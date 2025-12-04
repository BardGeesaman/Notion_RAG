"""
File parsing utilities for transcriptomics ingestion.

This module provides functions for parsing CSV/TSV files and extracting
gene identifiers with automatic column detection.
"""

from __future__ import annotations

from typing import Set, Tuple

import pandas as pd

from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_gene_set_from_file(file_path: str) -> Tuple[Set[str], pd.DataFrame]:
    """
    Load the DGE table and return a set of normalized gene identifiers and the full DataFrame.

    Args:
        file_path: Path to the transcriptomics file

    Returns:
        Tuple of (set of normalized genes, full DataFrame)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][TRANSCRIPTOMICS] File %s is empty", file_path)
        return set(), df

    # Detect gene identity column
    gene_column = None
    candidate_names = [
        "gene",
        "Gene",
        "gene_name",
        "Gene Symbol",
        "GeneSymbol",
        "symbol",
        "Symbol",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            gene_column = col
            break

    if not gene_column:
        # Try case-insensitive match
        for col in df.columns:
            col_lower = col.lower()
            if any(
                cand.lower() in col_lower or col_lower in cand.lower()
                for cand in candidate_names
            ):
                gene_column = col
                break

    if not gene_column:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Could not find gene identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find gene identity column in {file_path}"
        )

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Using column '%s' for gene identity in file %s",
        gene_column,
        file_path,
    )

    # Extract and normalize genes
    gene_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(gene_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_gene_identifier(raw_name_str)
        if normalized:
            gene_set.add(normalized)

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Extracted %d unique genes from %d rows in file %s",
        len(gene_set),
        total_rows,
        file_path,
    )

    return gene_set, df

