from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset, Feature
from scripts.dashboard.db_session import db_session


def _list_datasets(db: Session) -> List[Dict[str, str]]:
    datasets = (
        db.query(Dataset)
        .order_by(Dataset.updated_at.desc())
        .limit(200)  # avoid huge dropdowns
        .all()
    )
    return [
        {"id": str(d.id), "name": d.name, "omics_type": d.omics_type or "unknown"}
        for d in datasets
    ]


def _extract_fc_pval(feature: Feature) -> Optional[Dict[str, float]]:
    data = feature.external_ids or {}
    if not isinstance(data, dict):
        return None

    candidates = [data]
    for key in ("stats", "de_stats", "diffexp", "differential_expression"):
        nested = data.get(key)
        if isinstance(nested, dict):
            candidates.append(nested)

    def _first(keys):
        for obj in candidates:
            for k in keys:
                if k in obj and obj[k] is not None:
                    return obj[k]
        return None

    fc = _first(["log2FC", "log2fc", "log2_fold_change", "fold_change", "fc"])
    p = _first(["pvalue", "p_value", "adj_pvalue", "adj_p", "pval", "p_val"])

    try:
        fc = float(fc) if fc is not None else None
        p = float(p) if p is not None else None
    except Exception:
        return None

    if fc is None or p is None:
        return None
    return {"log2FC": fc, "pvalue": p}


def _get_volcano_data(db: Session, dataset_id: str, feature_type: str) -> pd.DataFrame:
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    if not dataset:
        return pd.DataFrame()

    rows = []
    for feat in dataset.features:
        if feature_type != "all" and (feat.feature_type or "").lower() != feature_type:
            continue
        stats = _extract_fc_pval(feat)
        if not stats:
            continue
        rows.append({"feature": feat.name, **stats})

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["-log10p"] = -np.log10(df["pvalue"])
    return df


def render() -> None:
    st.header("Volcano Plot")
    st.caption("Differential expression: log2 fold-change vs -log10(p-value)")

    with db_session() as db:
        datasets = _list_datasets(db)

    if not datasets:
        st.info("No datasets available in Postgres.")
        return

    dataset_labels = {f"{d['name']} ({d['omics_type']})": d["id"] for d in datasets}
    label = st.selectbox("Dataset", list(dataset_labels.keys()), index=0, key="volcano_dataset")
    dataset_id = dataset_labels[label]

    feature_type = st.selectbox(
        "Feature type",
        ["all", "gene", "protein", "metabolite", "lipid"],
        index=0,
        key="volcano_feature_type",
    )
    fc_thresh = st.slider("Fold-change threshold (abs)", 0.5, 3.0, 1.0, 0.1, key="volcano_fc_thresh")
    p_thresh = st.slider("p-value threshold", 0.0001, 0.1, 0.05, 0.0001, format="%.4f", key="volcano_p_thresh")

    if st.button("Refresh data", key="volcano_refresh"):
        st.rerun()

    if st.button("Generate Plot", key="volcano_generate"):
        with st.spinner("Computing volcano plot from Postgres..."):
            with db_session() as db:
                df = _get_volcano_data(db, dataset_id, feature_type)

        if df.empty:
            st.warning("No feature statistics found for this selection.")
            return

        df["significant"] = (df["pvalue"] <= p_thresh) & (np.abs(df["log2FC"]) >= fc_thresh)

        fig = px.scatter(
            df,
            x="log2FC",
            y="-log10p",
            color="significant",
            hover_name="feature",
            color_discrete_map={True: "#FF6B6B", False: "#4E79A7"},
            labels={"log2FC": "log2 Fold Change", "-log10p": "-log10(p-value)"},
        )
        fig.add_vline(x=fc_thresh, line_dash="dash", line_color="#999999")
        fig.add_vline(x=-fc_thresh, line_dash="dash", line_color="#999999")
        fig.add_hline(y=-np.log10(p_thresh), line_dash="dash", line_color="#999999")

        st.plotly_chart(fig, use_container_width=True)
        st.download_button(
            "Download CSV",
            data=df.to_csv(index=False),
            file_name="volcano_data.csv",
            mime="text/csv",
        )

