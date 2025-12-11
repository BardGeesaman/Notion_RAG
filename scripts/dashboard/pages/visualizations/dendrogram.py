from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.cluster import hierarchy as sch
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Feature
from scripts.dashboard.db_session import db_session


def _build_matrix(db: Session) -> pd.DataFrame:
    datasets: List[Dataset] = db.query(Dataset).all()
    if not datasets:
        return pd.DataFrame()

    feature_names = set()
    rows: Dict[str, Dict[str, float]] = {}

    for ds in datasets:
        row: Dict[str, float] = {}
        for feat in ds.features:
            # Use a simple presence/weight value if available
            val = 1.0
            if isinstance(feat, Feature) and feat.external_ids:
                # Try common keys for a numeric weight
                for key in ("log2fc", "pvalue", "weight"):
                    if key in feat.external_ids:
                        try:
                            val = float(feat.external_ids[key])
                            break
                        except Exception:
                            pass
            row[feat.name] = val
            feature_names.add(feat.name)
        rows[ds.name or str(ds.id)] = row

    if not feature_names:
        return pd.DataFrame()

    mat = pd.DataFrame(0.0, index=list(rows.keys()), columns=sorted(feature_names))
    for ds_label, row in rows.items():
        for fname, val in row.items():
            mat.at[ds_label, fname] = val
    return mat


def _make_dendrogram(mat: pd.DataFrame) -> go.Figure:
    # Compute distance and linkage
    data = mat.values
    if data.shape[0] <= 1:
        return go.Figure()

    # Use correlation distance
    corr = np.corrcoef(data)
    dist = 1 - corr
    # Ensure diagonals zero
    np.fill_diagonal(dist, 0)

    # Convert to condensed form for linkage
    condensed = sch.distance.squareform(dist, checks=False)
    Z = sch.linkage(condensed, method="average")

    dendro = sch.dendrogram(Z, labels=mat.index.tolist(), no_plot=True)
    icoord = np.array(dendro["icoord"])
    dcoord = np.array(dendro["dcoord"])
    labels = dendro["ivl"]

    fig = go.Figure()
    for xs, ys in zip(icoord, dcoord):
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color="#4E79A7", width=2)))

    fig.update_layout(
        title="Dataset Similarity Dendrogram",
        xaxis=dict(showticklabels=True, ticktext=labels, tickvals=range(5, 5 * len(labels) + 1, 10), tickangle=45),
        yaxis=dict(title="Distance"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=40, r=20, t=60, b=120),
    )
    return fig


def render() -> None:
    st.header("Dataset Similarity")
    st.caption("Hierarchical clustering of datasets based on shared features.")

    with db_session() as db:
        mat = _build_matrix(db)

    if mat.empty or mat.shape[0] <= 1:
        st.info("Not enough datasets or features to build a dendrogram.")
        return

    fig = _make_dendrogram(mat)
    st.plotly_chart(fig, use_container_width=True)

    st.download_button(
        "Download Matrix (CSV)",
        data=mat.to_csv(),
        file_name="dataset_similarity_matrix.csv",
        mime="text/csv",
    )


__all__ = ["render"]
