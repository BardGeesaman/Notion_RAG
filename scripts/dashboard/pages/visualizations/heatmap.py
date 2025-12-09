from __future__ import annotations

from typing import Dict, List

import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.cluster import hierarchy as sch
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Feature
from scripts.dashboard.db_session import db_session


def _list_datasets(db: Session) -> List[Dict[str, str]]:
    datasets = db.query(Dataset).order_by(Dataset.updated_at.desc()).limit(200).all()
    return [
        {"id": str(d.id), "name": d.name, "omics_type": d.omics_type or "unknown"}
        for d in datasets
    ]


def _feature_value(feature: Feature) -> float:
    ext = feature.external_ids or {}
    if isinstance(ext, dict):
        for key in ("log2FC", "log2fc", "value", "fold_change"):
            if key in ext and ext[key] is not None:
                try:
                    return float(ext[key])
                except Exception:
                    continue
    return 1.0  # presence indicator fallback


def _build_matrix(db: Session, dataset_ids: List[str]) -> pd.DataFrame:
    if not dataset_ids:
        return pd.DataFrame()

    datasets = db.query(Dataset).filter(Dataset.id.in_(dataset_ids)).all()
    if not datasets:
        return pd.DataFrame()

    feature_names = set()
    data: Dict[str, Dict[str, float]] = {}

    for ds in datasets:
        col = {}
        for feat in ds.features:
            val = _feature_value(feat)
            col[feat.name] = val
            feature_names.add(feat.name)
        data[str(ds.id)] = col

    if not feature_names:
        return pd.DataFrame()

    mat = pd.DataFrame(index=sorted(feature_names), columns=[str(ds.id) for ds in datasets])
    for ds in datasets:
        col = data.get(str(ds.id), {})
        mat[str(ds.id)] = mat.index.map(lambda f: col.get(f, 0.0))

    return mat


def _cluster_matrix(mat: pd.DataFrame) -> pd.DataFrame:
    if mat.shape[0] > 1:
        row_link = sch.linkage(mat.values, method="average")
        row_order = sch.leaves_list(row_link)
        mat = mat.iloc[row_order, :]
    if mat.shape[1] > 1:
        col_link = sch.linkage(mat.values.T, method="average")
        col_order = sch.leaves_list(col_link)
        mat = mat.iloc[:, col_order]
    return mat


def render() -> None:
    st.header("Heatmap")
    st.caption("Feature × Dataset intensity matrix.")

    with db_session() as db:
        datasets = _list_datasets(db)

    if not datasets:
        st.info("No datasets available.")
        return

    options = {f"{d['name']} ({d['omics_type']})": d["id"] for d in datasets}
    selected_labels = st.multiselect(
        "Datasets",
        list(options.keys()),
        default=list(options.keys())[:3],
        key="heatmap_datasets",
    )
    selected_ids = [options[lbl] for lbl in selected_labels]

    if st.button("Refresh data", key="heatmap_refresh"):
        st.rerun()

    if st.button("Generate Heatmap", key="heatmap_generate"):
        with st.spinner("Building heatmap from Postgres..."):
            with db_session() as db:
                mat = _build_matrix(db, selected_ids)

        if mat.empty:
            st.warning("No features found for selected datasets.")
            return

        clustered = _cluster_matrix(mat.fillna(0.0))
        fig = px.imshow(
            clustered,
            color_continuous_scale="Viridis",
            aspect="auto",
            labels=dict(color="Value"),
        )
        st.plotly_chart(fig, use_container_width=True)
        st.download_button(
            "Download CSV",
            data=clustered.to_csv(),
            file_name="heatmap_matrix.csv",
            mime="text/csv",
        )

