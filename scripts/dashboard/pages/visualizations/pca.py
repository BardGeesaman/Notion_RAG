from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Program
from scripts.dashboard.db_session import db_session


def _list_programs(db: Session) -> List[Dict[str, str]]:
    programs = db.query(Program).order_by(Program.updated_at.desc()).limit(100).all()
    return [{"id": str(p.id), "name": p.name} for p in programs]


def _list_datasets(db: Session) -> List[Dict[str, str]]:
    datasets = db.query(Dataset).order_by(Dataset.updated_at.desc()).limit(200).all()
    return [
        {"id": str(d.id), "name": d.name, "omics_type": d.omics_type or "unknown"}
        for d in datasets
    ]


def _load_dataset_info(db: Session, dataset_ids: List[str]) -> List[Dict[str, Any]]:
    """Materialize datasets and feature vectors inside the session."""
    if not dataset_ids:
        return []
    objs = db.query(Dataset).filter(Dataset.id.in_(dataset_ids)).limit(200).all()
    materialized: List[Dict[str, Any]] = []
    for ds in objs:
        feat_map: Dict[str, float] = {}
        for feat in ds.features:
            val = 1.0
            ext = feat.external_ids or {}
            if isinstance(ext, dict):
                for key in ("log2FC", "log2fc", "value", "fold_change"):
                    if key in ext and ext[key] is not None:
                        try:
                            val = float(ext[key])
                        except Exception:
                            pass
                        break
            feat_map[feat.name] = val
        materialized.append(
            {
                "id": str(ds.id),
                "name": ds.name,
                "omics_type": ds.omics_type or "unknown",
                "disease": ds.disease or [],
                "programs": [p.name for p in ds.programs] if ds.programs else [],
                "features": feat_map,
            }
        )
    return materialized


def _build_matrix(datasets: List[Dict[str, Any]]) -> pd.DataFrame:
    if not datasets:
        return pd.DataFrame()
    feature_names = set()
    for ds in datasets:
        feature_names.update(ds.get("features", {}).keys())
    if not feature_names:
        return pd.DataFrame()

    matrix = pd.DataFrame(index=[ds["id"] for ds in datasets], columns=sorted(feature_names))
    for ds in datasets:
        vec = ds.get("features", {})
        matrix.loc[ds["id"]] = [vec.get(f, 0.0) for f in matrix.columns]
    return matrix.astype(float)


def render() -> None:
    st.header("PCA Scatter")
    st.caption("Dataset similarity via PCA (Postgres-backed).")

    with db_session() as db:
        programs = _list_programs(db)
        datasets = _list_datasets(db)

    color_by = st.selectbox("Color by", ["program", "omics_type", "disease"], index=2, key="pca_color_by")

    program_options = {"(none)": None}
    program_options.update({p["name"]: p["id"] for p in programs})
    program_choice = st.selectbox("Program (optional)", list(program_options.keys()), index=0, key="pca_program")
    program_id = program_options[program_choice]

    dataset_options = {f"{d['name']} ({d['omics_type']})": d["id"] for d in datasets}

    if program_id:
        # filter datasets in program inside session and materialize IDs
        with db_session() as db:
            program = db.query(Program).filter(Program.id == program_id).first()
            ds_filtered = program.datasets if program else []
            dataset_options = {f"{d.name} ({d.omics_type})": str(d.id) for d in ds_filtered}

    selected_labels = st.multiselect(
        "Datasets",
        list(dataset_options.keys()),
        default=list(dataset_options.keys())[:5],
        key="pca_datasets",
    )
    selected_ids = {dataset_options[lbl] for lbl in selected_labels}

    if st.button("Refresh data", key="pca_refresh"):
        st.rerun()

    if st.button("Run PCA", key="pca_generate"):
        with st.spinner("Running PCA from Postgres data..."):
            with db_session() as db:
                ds_info = _load_dataset_info(db, list(selected_ids))
                mat = _build_matrix(ds_info)

        if mat.empty:
            st.warning("No feature data found for selected datasets.")
            return

        pca = PCA(n_components=min(2, mat.shape[1]))
        coords = pca.fit_transform(mat)
        df = pd.DataFrame(coords, columns=["PC1", "PC2"][: coords.shape[1]])
        df["dataset_id"] = mat.index

        # color label using materialized info
        info_lookup = {d["id"]: d for d in ds_info}
        labels = []
        for ds_id in df["dataset_id"]:
            ds = info_lookup.get(ds_id)
            if not ds:
                labels.append("unknown")
            elif color_by == "omics_type":
                labels.append(ds.get("omics_type") or "unknown")
            elif color_by == "disease":
                labels.append(", ".join(ds.get("disease") or []) or "unknown")
            else:
                programs = ds.get("programs") or []
                labels.append(programs[0] if programs else "unknown")
        df["color"] = labels

        fig = px.scatter(
            df,
            x="PC1",
            y="PC2" if "PC2" in df.columns else "PC1",
            hover_name="dataset_id",
            color="color",
            labels={
                "PC1": f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)",
                "PC2": f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)" if pca.n_components_ > 1 else "PC2",
                "color": color_by,
            },
        )
        st.plotly_chart(fig, width='stretch')
        st.download_button(
            "Download PCA CSV",
            data=df.to_csv(index=False),
            file_name="pca_coordinates.csv",
            mime="text/csv",
        )

