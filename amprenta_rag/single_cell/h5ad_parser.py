"""h5ad loader + basic validation for single-cell ingestion."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def load_h5ad(path: str):
    """Load an .h5ad file into an AnnData object."""
    try:
        import anndata as ad  # type: ignore[import-not-found]
    except Exception as e:  # noqa: BLE001
        raise ImportError('load_h5ad requires "anndata" (and h5py)') from e

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    return ad.read_h5ad(str(p))


def extract_metadata(adata) -> Dict[str, Any]:
    """Extract basic metadata from AnnData."""
    n_cells = int(getattr(adata, "n_obs", 0))
    n_genes = int(getattr(adata, "n_vars", 0))

    var_names = list(getattr(getattr(adata, "var_names", []), "tolist", lambda: list(adata.var_names))())
    obs_names = list(getattr(getattr(adata, "obs_names", []), "tolist", lambda: list(adata.obs_names))())

    return {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "var_names": var_names,
        "obs_names": obs_names,
    }


def validate_h5ad(adata) -> None:
    """Validate minimal AnnData structure required by pipeline."""
    if adata is None:
        raise ValueError("AnnData is None")
    if not hasattr(adata, "X"):
        raise ValueError("AnnData missing X matrix")
    if getattr(adata, "n_obs", 0) <= 0:
        raise ValueError("AnnData has no observations (cells)")
    if getattr(adata, "n_vars", 0) <= 0:
        raise ValueError("AnnData has no variables (genes)")
    if not hasattr(adata, "obs") or not hasattr(adata, "var"):
        raise ValueError("AnnData missing obs/var")


__all__ = ["load_h5ad", "extract_metadata", "validate_h5ad"]


