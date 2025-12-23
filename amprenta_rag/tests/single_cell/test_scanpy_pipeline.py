from __future__ import annotations

import numpy as np
import pytest

from amprenta_rag.single_cell.scanpy_pipeline import run_preprocessing


def test_run_preprocessing_synthetic_anndata():
    pytest.importorskip("anndata")
    pytest.importorskip("scanpy")

    import anndata as ad  # type: ignore[import-not-found]

    rng = np.random.default_rng(0)
    X = rng.poisson(lam=1.0, size=(100, 200)).astype(np.float32)
    adata = ad.AnnData(X)
    adata.var_names = [f"GENE{i}" for i in range(200)]
    adata.obs_names = [f"CELL{i}" for i in range(100)]

    out = run_preprocessing(adata)
    assert "X_umap" in out.obsm
    assert "leiden" in out.obs


def test_resolution_affects_clustering():
    pytest.importorskip("anndata")
    pytest.importorskip("scanpy")

    import anndata as ad  # type: ignore[import-not-found]

    rng = np.random.default_rng(42)
    n_cells = 180
    n_genes = 300

    # Create 3 clearly-separated groups in gene space.
    labels = np.repeat([0, 1, 2], n_cells // 3)
    X = rng.uniform(0.1, 2.0, size=(n_cells, n_genes)).astype(np.float32)
    X[labels == 0, :50] += 5.0
    X[labels == 1, 50:100] += 5.0
    X[labels == 2, 100:150] += 5.0

    adata = ad.AnnData(X)
    adata.var_names = [f"GENE{i}" for i in range(n_genes)]
    adata.obs_names = [f"CELL{i}" for i in range(n_cells)]

    low = run_preprocessing(adata, resolution=0.1)
    high = run_preprocessing(adata, resolution=1.0)

    low_k = int(low.obs["leiden"].nunique())
    high_k = int(high.obs["leiden"].nunique())
    # Higher resolution should not reduce the number of clusters.
    assert high_k >= low_k
    assert high_k != low_k


