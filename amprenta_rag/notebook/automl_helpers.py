"""Helpers for AutoML notebook templates.

These helpers are designed for use inside JupyterHub notebooks and for headless
Papermill execution.
"""

from __future__ import annotations
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
from uuid import UUID

import httpx
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)

from amprenta_rag.ml.registry import get_registry


API_BASE = os.environ.get("API_URL", "http://localhost:8000").rstrip("/")


def _get_dataset_metadata(dataset_id: str) -> dict:
    with httpx.Client(timeout=60) as client:
        r = client.get(f"{API_BASE}/api/v1/datasets/{dataset_id}")
    r.raise_for_status()
    return dict(r.json() or {})


def _read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))
    suf = p.suffix.lower()
    if suf in {".csv"}:
        return pd.read_csv(p)
    if suf in {".tsv", ".txt"}:
        return pd.read_csv(p, sep="\t")
    if suf in {".parquet"}:
        return pd.read_parquet(p)
    if suf in {".xlsx"}:
        return pd.read_excel(p)
    raise ValueError(f"Unsupported dataset file type: {suf}")


def load_dataset_as_dataframe(dataset_id: str) -> pd.DataFrame:
    """Load a dataset into a DataFrame.

    Strategy:
    - Fetch dataset metadata from API (/api/v1/datasets/{id})
    - Read first path from dataset.file_paths on local disk (JupyterHub shared volume)
    - Fallback: if dataset.file_urls exists and is http(s), attempt to fetch CSV
    """
    meta = _get_dataset_metadata(str(dataset_id))
    file_paths = list(meta.get("file_paths") or [])
    if file_paths:
        return _read_table(str(file_paths[0]))

    file_urls = list(meta.get("file_urls") or [])
    if file_urls:
        u = str(file_urls[0])
        if u.startswith("http://") or u.startswith("https://"):
            with httpx.Client(timeout=120) as client:
                r = client.get(u)
            r.raise_for_status()
            # best-effort CSV
            return pd.read_csv(pd.io.common.BytesIO(r.content))

    raise FileNotFoundError("Dataset has no accessible file_paths/file_urls")


def register_trained_model(
    model: Any,
    name: str,
    metrics: Dict[str, float],
    dataset_id: str,
    *,
    version: str = "1.0.0",
    model_type: str = "automl",
    framework: str = "sklearn",
    features: Optional[List[str]] = None,
    hyperparameters: Optional[Dict[str, Any]] = None,
    description: Optional[str] = None,
) -> Any:
    """Register a trained model in MLModelRegistry (artifact + DB row)."""
    reg = get_registry()
    return reg.register_model(
        name=name,
        version=version,
        model_type=model_type,
        framework=framework,
        model_object=model,
        features=features,
        hyperparameters=hyperparameters,
        metrics=metrics,
        description=description,
        training_dataset_id=UUID(str(dataset_id)),
    )


def generate_classification_report(y_true: Sequence[float], y_pred: Sequence[float], y_proba: Optional[Sequence[float]] = None) -> Dict[str, float]:
    yt = np.asarray(y_true)
    yp = np.asarray(y_pred)
    out: Dict[str, float] = {"accuracy": float(accuracy_score(yt, yp))}
    if y_proba is not None:
        out["auc"] = float(roc_auc_score(yt, np.asarray(y_proba)))
    return out


def generate_regression_report(y_true: Sequence[float], y_pred: Sequence[float]) -> Dict[str, float]:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    rmse = float(np.sqrt(mean_squared_error(yt, yp)))
    return {
        "rmse": rmse,
        "mae": float(mean_absolute_error(yt, yp)),
        "r2": float(r2_score(yt, yp)),
    }


def plot_confusion_matrix(y_true: Sequence[float], y_pred: Sequence[float]):  # noqa: ANN001
    import matplotlib.pyplot as plt

    cm = confusion_matrix(np.asarray(y_true), np.asarray(y_pred))
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(cm, cmap="Blues")
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center")
    return fig


def plot_elbow_curve(X, max_k: int = 10):  # noqa: ANN001
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans

    inertias = []
    ks = list(range(1, int(max_k) + 1))
    for k in ks:
        km = KMeans(n_clusters=k, n_init="auto", random_state=42)
        km.fit(X)
        inertias.append(float(km.inertia_))

    fig, ax = plt.subplots()
    ax.plot(ks, inertias, marker="o")
    ax.set_xlabel("k")
    ax.set_ylabel("Inertia")
    ax.set_title("Elbow curve")
    return fig


__all__ = [
    "load_dataset_as_dataframe",
    "register_trained_model",
    "generate_classification_report",
    "generate_regression_report",
    "plot_confusion_matrix",
    "plot_elbow_curve",
]


