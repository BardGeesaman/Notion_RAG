"""
File parsing utilities for lipidomics ingestion.

This module provides functions for parsing CSV/TSV files and extracting
lipid species with automatic column detection.
"""

from __future__ import annotations

from typing import Set, Tuple

import pandas as pd

from amprenta_rag.ingestion.lipidomics.normalization import normalize_lipid_species
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_species_from_file(file_path: str) -> Tuple[Set[str], int]:
    """
    Extract and normalize lipid species from a CSV/TSV file.

    Args:
        file_path: Path to the lipidomics file

    Returns:
        Tuple of (set of normalized species, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][LIPIDOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][LIPIDOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][LIPIDOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect lipid identity column
    lipid_column = None
    candidate_names = ["species", "lipid", "Lipid", "Name", "Molecule", "compound", "metabolite"]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            lipid_column = col
            break

    if not lipid_column:
        # Try case-insensitive match
        for col in df.columns:
            if any(cand.lower() in col.lower() for cand in candidate_names):
                lipid_column = col
                break

    if not lipid_column:
        logger.error(
            "[INGEST][LIPIDOMICS] Could not find lipid identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(f"Could not find lipid identity column in {file_path}")

    logger.info(
        "[INGEST][LIPIDOMICS] Using column '%s' for lipid identity in file %s",
        lipid_column,
        file_path,
    )

    # Extract and normalize species
    species_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(lipid_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_lipid_species(raw_name_str)
        if normalized:
            species_set.add(normalized)
        else:
            # Keep raw name if normalization fails (for now)
            species_set.add(raw_name_str)

    logger.info(
        "[INGEST][LIPIDOMICS] Extracted %d unique species from %d rows in file %s",
        len(species_set),
        total_rows,
        file_path,
    )

    return species_set, total_rows

