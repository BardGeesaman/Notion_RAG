"""Scanpy preprocessing pipeline for single-cell AnnData objects."""

from __future__ import annotations


def run_preprocessing(adata, resolution: float = 0.5):
    """Run preprocessing on an AnnData object.

    Pipeline:
    - QC filtering: filter_cells(min_genes=200), filter_genes(min_cells=3)
    - Normalization: normalize_total(1e4), log1p()
    - HVG: highly_variable_genes(n_top_genes=2000)
    - PCA: pca(n_comps=50)
    - Neighbors: neighbors(n_neighbors=10, n_pcs=40)
    - UMAP: umap()
    - Clustering: leiden(resolution=0.5)

    Returns AnnData with:
    - adata.obsm['X_umap']
    - adata.obs['leiden']
    """
    try:
        import scanpy as sc  # type: ignore[import-not-found]
    except Exception as e:  # noqa: BLE001
        raise ImportError('run_preprocessing requires "scanpy"') from e

    # Work on a copy to avoid mutating caller unexpectedly.
    ad = adata.copy()

    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    n_top = min(2000, int(getattr(ad, "n_vars", 0) or 0))
    if n_top <= 0:
        raise ValueError("AnnData has no genes after QC filtering")
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top, subset=True)
    sc.pp.scale(ad, max_value=10)

    n_vars = int(getattr(ad, "n_vars", 0) or 0)
    n_obs = int(getattr(ad, "n_obs", 0) or 0)
    n_comps = min(50, max(2, min(n_vars, max(2, n_obs - 1))))
    sc.tl.pca(ad, n_comps=n_comps, svd_solver="arpack")
    n_pcs = min(40, n_comps)
    sc.pp.neighbors(ad, n_neighbors=10, n_pcs=n_pcs)
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=float(resolution))

    return ad


__all__ = ["run_preprocessing"]


