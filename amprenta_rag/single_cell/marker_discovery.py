"""Marker discovery utilities for single-cell clusters."""

from __future__ import annotations

import pandas as pd


def find_markers(adata, groupby: str = "leiden", *, top_n: int = 50) -> pd.DataFrame:
    """Find top markers per cluster using Scanpy rank_genes_groups (wilcoxon).

    Returns a DataFrame with columns:
    - gene_symbol, cluster_id, log2_fold_change, pval, pval_adj, pct_in_cluster, pct_out_cluster
    """
    try:
        import scanpy as sc  # type: ignore[import-not-found]
        import numpy as np  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError('find_markers requires "scanpy"') from e

    if groupby not in adata.obs:
        raise ValueError(f"Missing adata.obs['{groupby}']")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon")
    df = sc.get.rank_genes_groups_df(adata, group=None)

    # Normalize column names
    rename = {
        "names": "gene_symbol",
        "group": "cluster_id",
        "logfoldchanges": "log2_fold_change",
        "pvals": "pval",
        "pvals_adj": "pval_adj",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Compute pct_in/out: fraction of cells with expression > 0 for each gene in/out cluster.
    # This is a best-effort approximation for MVP.
    # Ensure var_names are strings for indexing
    var_names = [str(x) for x in adata.var_names]

    def _pct_for_gene(gene: str, mask) -> float:
        if gene not in var_names:
            return float("nan")
        idx = var_names.index(gene)
        x = adata.X[:, idx]
        try:
            # handle sparse
            arr = x.toarray().ravel()
        except Exception:
            arr = np.asarray(x).ravel()
        sel = arr[mask]
        if sel.size == 0:
            return float("nan")
        return float((sel > 0).mean())

    pct_in = []
    pct_out = []
    groups = adata.obs[groupby].astype(str)
    for _, row in df.iterrows():
        gene = str(row["gene_symbol"])
        cluster = str(row["cluster_id"])
        in_mask = (groups == cluster).to_numpy()
        out_mask = (groups != cluster).to_numpy()
        pct_in.append(_pct_for_gene(gene, in_mask))
        pct_out.append(_pct_for_gene(gene, out_mask))

    df["pct_in_cluster"] = pct_in
    df["pct_out_cluster"] = pct_out

    # Top N per cluster
    if "cluster_id" in df.columns:
        df = df.sort_values(["cluster_id", "pval_adj"], ascending=[True, True])
        df = df.groupby("cluster_id", as_index=False).head(int(top_n))

    keep = [
        "gene_symbol",
        "cluster_id",
        "log2_fold_change",
        "pval",
        "pval_adj",
        "pct_in_cluster",
        "pct_out_cluster",
    ]
    return df[[c for c in keep if c in df.columns]].reset_index(drop=True)


__all__ = ["find_markers"]


