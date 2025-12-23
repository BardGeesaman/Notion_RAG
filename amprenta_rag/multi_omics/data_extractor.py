"""Multi-omics data extraction helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List
from uuid import UUID

import pandas as pd

from amprenta_rag.database.models import Dataset as DatasetModel


def _try_read_matrix(file_path: str) -> pd.DataFrame:
    """Read a tabular omics matrix from disk.

    Expected: first column is feature IDs (or a known feature column), remaining columns are samples.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"Dataset file not found: {file_path}")

    for sep in (None, "\t", ",", ";"):
        try:
            df = pd.read_csv(p, sep=sep, engine="python")
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    raise ValueError(f"Unable to parse dataset file as a matrix: {file_path}")


def _detect_feature_column(df: pd.DataFrame) -> str:
    for cand in ("feature", "feature_id", "gene", "gene_id", "id", "name"):
        if cand in df.columns:
            return cand
    return df.columns[0]


def _as_feature_by_sample(df: pd.DataFrame) -> pd.DataFrame:
    feat_col = _detect_feature_column(df)
    df2 = df.copy().set_index(feat_col)
    num = df2.select_dtypes(include=["number"])
    if num.empty:
        num = df2.apply(pd.to_numeric, errors="coerce")
        num = num.dropna(axis=1, how="all")
    if num.empty:
        raise ValueError("No numeric sample columns detected in matrix")
    return num


def extract_feature_matrix(dataset_id: UUID, db) -> pd.DataFrame:
    """Extract a features x samples numeric matrix from Dataset.file_paths[0]."""
    ds: DatasetModel | None = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    if not ds:
        raise ValueError("Dataset not found")
    fp = (ds.file_paths or [None])[0]
    if not fp:
        raise ValueError(f"Dataset {ds.id} has no file_paths")
    raw = _try_read_matrix(str(fp))
    return _as_feature_by_sample(raw)


def align_samples(matrices: Dict[str, pd.DataFrame], sample_mapping: List[dict]) -> Dict[str, pd.DataFrame]:
    """Align sample columns across omics matrices using a provided mapping.

    Supported `sample_mapping` shapes:

    A) Row-oriented (recommended):
      [
        {"sample_id": "S1", "transcriptomics": "RNA_A", "proteomics": "PROT_A"},
        {"sample_id": "S2", "transcriptomics": "RNA_B", "proteomics": "PROT_B"},
      ]
      - For each omics_type key present in `matrices`, the value must match a column in that matrix.
      - Output matrices are reordered and columns renamed to unified `sample_id`.

    B) Long-form:
      [
        {"sample_id": "S1", "omics_type": "transcriptomics", "source_sample_id": "RNA_A"},
        {"sample_id": "S1", "omics_type": "proteomics", "source_sample_id": "PROT_A"},
        ...
      ]
    """
    if not sample_mapping:
        return dict(matrices)

    # Determine if row-oriented mapping is present
    has_row_oriented = any(isinstance(m, dict) and "sample_id" in m for m in sample_mapping) and any(
        any(k in (m or {}) for k in matrices.keys()) for m in sample_mapping
    )

    if has_row_oriented:
        order = [str(m.get("sample_id")) for m in sample_mapping if m.get("sample_id")]
        out: Dict[str, pd.DataFrame] = {}
        for omics_type, df in matrices.items():
            src_cols: List[str] = []
            for m in sample_mapping:
                if m.get("sample_id") is None:
                    continue
                src = m.get(omics_type)
                if src is None:
                    raise ValueError(f"sample_mapping missing '{omics_type}' for sample_id={m.get('sample_id')}")
                src_cols.append(str(src))
            missing = [c for c in src_cols if c not in df.columns]
            if missing:
                raise ValueError(f"{omics_type} matrix missing mapped sample columns: {missing[:5]}")
            aligned = df.loc[:, src_cols].copy()
            aligned.columns = order
            out[omics_type] = aligned
        return out

    # Long-form mapping
    order: List[str] = []
    per_omics: Dict[str, Dict[str, str]] = {k: {} for k in matrices.keys()}
    for m in sample_mapping:
        if not isinstance(m, dict):
            continue
        sid = m.get("sample_id")
        if sid and sid not in order:
            order.append(str(sid))
        ot = m.get("omics_type")
        src = m.get("source_sample_id") or m.get("source") or m.get("sample")
        if ot in per_omics and sid and src:
            per_omics[str(ot)][str(sid)] = str(src)

    out2: Dict[str, pd.DataFrame] = {}
    for omics_type, df in matrices.items():
        mapping = per_omics.get(omics_type) or {}
        src_cols = []
        for sid in order:
            if sid not in mapping:
                raise ValueError(f"sample_mapping missing mapping for omics_type={omics_type}, sample_id={sid}")
            src_cols.append(mapping[sid])
        missing = [c for c in src_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{omics_type} matrix missing mapped sample columns: {missing[:5]}")
        aligned = df.loc[:, src_cols].copy()
        aligned.columns = order
        out2[omics_type] = aligned

    return out2


__all__ = ["extract_feature_matrix", "align_samples"]


