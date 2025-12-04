"""
File parsing utilities for metabolomics ingestion.

This module provides functions for parsing CSV/TSV files and extracting
metabolite names with automatic column detection.
"""

from __future__ import annotations

from typing import Set, Tuple

import pandas as pd

from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_metabolite_set_from_file(file_path: str) -> Tuple[Set[str], int]:
    """
    Load CSV/TSV and return a set of normalized metabolite names.

    Args:
        file_path: Path to the metabolomics file

    Returns:
        Tuple of (set of normalized metabolites, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][METABOLOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][METABOLOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][METABOLOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect metabolite identity column
    metabolite_column = None
    candidate_names = [
        "metabolite",
        "Metabolite",
        "compound",
        "Compound",
        "name",
        "Name",
        "molecule",
        "Molecule",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            metabolite_column = col
            break

    if not metabolite_column:
        # Try case-insensitive match
        for col in df.columns:
            if any(cand.lower() in col.lower() for cand in candidate_names):
                metabolite_column = col
                break

    if not metabolite_column:
        logger.error(
            "[INGEST][METABOLOMICS] Could not find metabolite identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find metabolite identity column in {file_path}"
        )

    logger.info(
        "[INGEST][METABOLOMICS] Using column '%s' for metabolite identity in file %s",
        metabolite_column,
        file_path,
    )

    # Extract and normalize metabolites
    metabolite_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(metabolite_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_metabolite_name(raw_name_str)
        if normalized:
            metabolite_set.add(normalized)

    logger.info(
        "[INGEST][METABOLOMICS] Extracted %d unique metabolites from %d rows in file %s",
        len(metabolite_set),
        total_rows,
        file_path,
    )

    return metabolite_set, total_rows

