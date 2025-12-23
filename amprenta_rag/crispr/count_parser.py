"""CRISPR count matrix parsing utilities."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def parse_count_matrix(path: str) -> pd.DataFrame:
    """Parse a MAGeCK-style count matrix.

    Expected columns:
    - sgRNA
    - Gene
    - one or more sample columns
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    # MAGeCK count tables are commonly TSV, but allow CSV as well.
    # Use sep=None to infer delimiter via Python engine.
    df = pd.read_csv(p, sep=None, engine="python")
    validate_count_matrix(df)
    return df


def validate_count_matrix(df: pd.DataFrame) -> None:
    """Validate required columns for a MAGeCK count matrix."""
    cols = {c.strip(): c for c in df.columns}
    if "sgRNA" not in cols:
        raise ValueError("Count matrix missing required column: sgRNA")
    if "Gene" not in cols:
        raise ValueError("Count matrix missing required column: Gene")
    if len(df.columns) < 3:
        raise ValueError("Count matrix must include at least one sample column in addition to sgRNA and Gene")


__all__ = ["parse_count_matrix", "validate_count_matrix"]


