"""Batch effect correction utilities (ComBat MVP)."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

import pandas as pd

from amprenta_rag.database.models import BatchCorrection as BatchCorrectionModel
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _try_read_matrix(file_path: str) -> pd.DataFrame:
    """Read a tabular omics matrix from disk.

    Expected: first column is feature IDs, remaining columns are numeric samples.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"Dataset file not found: {file_path}")

    # Try common separators.
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


def _as_numeric_matrix(df: pd.DataFrame, dataset_prefix: str) -> pd.DataFrame:
    feat_col = _detect_feature_column(df)
    df2 = df.copy()
    df2 = df2.set_index(feat_col)
    # keep numeric columns
    num = df2.select_dtypes(include=["number"])
    if num.empty:
        # attempt to coerce
        num = df2.apply(pd.to_numeric, errors="coerce")
        num = num.dropna(axis=1, how="all")
    if num.empty:
        raise ValueError("No numeric sample columns detected in matrix")
    # prefix columns to avoid collisions across datasets
    num.columns = [f"{dataset_prefix}:{c}" for c in num.columns]
    return num


def _batch_labels_for_columns(
    columns: List[str], batch_map: Dict[str, Any]
) -> List[str]:
    """Resolve a batch label per sample column.

    Supports:
    - exact sample mapping: {"ds:sampleA": "batch1"}
    - dataset-wide mapping: {"<dataset_uuid>": "batch1"} where columns start with "<dataset_uuid>:"
    """
    labels: List[str] = []
    for col in columns:
        if col in batch_map:
            labels.append(str(batch_map[col]))
            continue
        ds_prefix = col.split(":", 1)[0]
        if ds_prefix in batch_map:
            labels.append(str(batch_map[ds_prefix]))
            continue
        labels.append("unknown")
    return labels


def _combat_correct(df: pd.DataFrame, batch: List[str]) -> pd.DataFrame:
    """Run ComBat via optional dependency."""
    # Several python packages exist; try a few import paths.
    try:
        from combat.pycombat import pycombat  # type: ignore

        corrected = pycombat(df, batch=batch)
        return corrected
    except Exception:
        pass
    try:
        from pycombat.pycombat import pycombat  # type: ignore

        corrected = pycombat(df, batch=batch)
        return corrected
    except Exception:
        pass
    try:
        from pycombat import pycombat  # type: ignore

        corrected = pycombat(df, batch=batch)
        return corrected
    except Exception:
        pass
    raise ImportError("ComBat dependency not installed. Install combat>=0.3 or pycombat.")


def correct_batch_effects_df(
    matrix: pd.DataFrame,
    batch_labels: List[str],
    method: str = "combat",
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Correct batch effects for a provided matrix (features x samples)."""
    if matrix.empty:
        raise ValueError("matrix is empty")
    if len(batch_labels) != matrix.shape[1]:
        raise ValueError("batch_labels length must match number of sample columns")

    method_l = (method or "combat").lower()
    before = matrix.copy()

    if method_l == "combat":
        corrected = _combat_correct(matrix, batch=batch_labels)
    else:
        raise ValueError(f"Unknown method: {method}")

    stats = {
        "method": method_l,
        "n_features": int(matrix.shape[0]),
        "n_samples": int(matrix.shape[1]),
        "batches": {b: int(batch_labels.count(b)) for b in sorted(set(batch_labels))},
        "mean_before": float(before.to_numpy().mean()),
        "mean_after": float(corrected.to_numpy().mean()),
    }
    return corrected, stats


def correct_batch_effects(
    datasets: List[UUID],
    batch_map: Dict[str, Any],
    method: str = "combat",
) -> Dict[str, Any]:
    """Correct batch effects across datasets using their stored file paths.

    Returns:
      {
        "corrected_df": pd.DataFrame,
        "stats": dict,
        "corrected_dataset_id": UUID | None
      }
    """
    if not datasets:
        raise ValueError("datasets must be non-empty")

    # Load dataset records and matrices
    mats: List[pd.DataFrame] = []
    ds_names: List[str] = []
    with db_session() as db:
        ds_rows: List[DatasetModel] = (
            db.query(DatasetModel).filter(DatasetModel.id.in_(datasets)).all()
        )
        if len(ds_rows) != len(set(datasets)):
            raise ValueError("One or more datasets not found")
        for ds in ds_rows:
            fp = (ds.file_paths or [None])[0]
            if not fp:
                raise ValueError(f"Dataset {ds.id} has no file_paths")
            raw = _try_read_matrix(fp)
            mat = _as_numeric_matrix(raw, dataset_prefix=str(ds.id))
            mats.append(mat)
            ds_names.append(ds.name or str(ds.id))

    # Merge features by index; outer join, then fill missing with 0
    combined = pd.concat(mats, axis=1, join="outer").fillna(0.0)
    batch_labels = _batch_labels_for_columns(list(combined.columns), batch_map=batch_map)

    corrected, stats = correct_batch_effects_df(combined, batch_labels=batch_labels, method=method)

    corrected_dataset_id: Optional[UUID] = None
    try:
        out_dir = Path("data/batch_corrected")
        out_dir.mkdir(parents=True, exist_ok=True)
        corrected_path = out_dir / f"batch_corrected_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.csv"
        corrected.to_csv(corrected_path)

        with db_session() as db:
            # Create a new dataset entry pointing to corrected file (MVP).
            first = db.query(DatasetModel).filter(DatasetModel.id == datasets[0]).first()
            omics_type = getattr(first, "omics_type", "transcriptomics") if first else "transcriptomics"

            new_ds = DatasetModel(
                name=f"Batch-corrected ({method})",
                omics_type=omics_type,
                description=f"Batch-corrected from {len(datasets)} dataset(s).",
                file_paths=[str(corrected_path)],
                ingestion_status="complete",
            )
            db.add(new_ds)
            db.commit()
            db.refresh(new_ds)
            corrected_dataset_id = new_ds.id

            rec = BatchCorrectionModel(
                method=str(method),
                batch_map=batch_map,
                corrected_dataset_id=corrected_dataset_id,
            )
            db.add(rec)
            db.commit()
    except Exception as e:  # noqa: BLE001
        logger.warning("Failed to persist batch correction result: %s", e)

    return {"corrected_df": corrected, "stats": stats, "corrected_dataset_id": corrected_dataset_id}


__all__ = ["correct_batch_effects", "correct_batch_effects_df"]


