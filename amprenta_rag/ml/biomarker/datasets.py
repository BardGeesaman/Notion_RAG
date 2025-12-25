"""Dataset loading utilities for biomarker discovery."""

from __future__ import annotations

from typing import List, Tuple
from uuid import UUID

import numpy as np

from amprenta_rag.database.models import Experiment as ExperimentModel
from amprenta_rag.database.session import db_session
from amprenta_rag.multi_omics.data_extractor import extract_feature_matrix


def load_biomarker_dataset(
    experiment_id: UUID | str,
    group1_samples: List[str],
    group2_samples: List[str],
    omics_type: str,
    min_coverage: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Load a biomarker dataset for an experiment and two sample groups.

    Returns:
        X: (n_samples, n_features) float array
        y: (n_samples,) int array (0=group1, 1=group2)
        feature_names: list[str] length n_features
    """
    exp_id = UUID(str(experiment_id)) if not isinstance(experiment_id, UUID) else experiment_id
    g1 = [str(s) for s in group1_samples or []]
    g2 = [str(s) for s in group2_samples or []]
    if not g1 or not g2:
        raise ValueError("group1_samples and group2_samples must be non-empty")

    omics = str(omics_type or "").strip()
    if not omics:
        raise ValueError("omics_type is required")

    with db_session() as db:
        exp = db.query(ExperimentModel).filter(ExperimentModel.id == exp_id).first()
        if not exp:
            raise ValueError("Experiment not found")

        # Select first dataset matching requested omics_type
        ds_match = None
        for ds in getattr(exp, "datasets", []) or []:
            if str(getattr(ds, "omics_type", "") or "").lower() == omics.lower():
                ds_match = ds
                break
        if not ds_match:
            raise ValueError(f"No dataset found for omics_type='{omics}' on experiment {exp_id}")

        # features x samples
        mat = extract_feature_matrix(ds_match.id, db)

    # Align/limit to requested samples (column names are sample IDs in file)
    ordered_samples = [s for s in (g1 + g2) if s in mat.columns]
    if not ordered_samples:
        raise ValueError("None of the provided samples exist in the dataset matrix")

    mat2 = mat.loc[:, ordered_samples].copy()

    # Drop features with insufficient coverage
    cov = mat2.notna().mean(axis=1)
    mat2 = mat2.loc[cov >= float(min_coverage)]

    # Median impute per feature
    med = mat2.median(axis=1, skipna=True)
    mat2 = mat2.T  # samples x features
    mat2 = mat2.fillna(med)

    X = mat2.to_numpy(dtype=float)
    feature_names = [str(x) for x in mat2.columns.tolist()]

    y = np.array([0] * len([s for s in g1 if s in ordered_samples]) + [1] * len([s for s in g2 if s in ordered_samples]), dtype=int)
    if X.shape[0] != y.shape[0]:
        # This can happen if some group sample IDs were missing; rebuild y from ordered_samples
        y = np.array([0 if s in g1 else 1 for s in ordered_samples], dtype=int)

    return X, y, feature_names


__all__ = ["load_biomarker_dataset"]


