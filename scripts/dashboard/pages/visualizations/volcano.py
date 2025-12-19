from __future__ import annotations


import numpy as np
import pandas as pd
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Dataset
from scripts.dashboard.db_session import db_session
from scripts.dashboard.utils.viz_helpers import (
    create_volcano_plot,
    extract_feature_stats,
    list_datasets_for_dropdown,
)


def _get_volcano_data(db: Session, dataset_id: str, feature_type: str) -> pd.DataFrame:
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    if not dataset:
        return pd.DataFrame()

    rows = []
    for feat in dataset.features:
        if feature_type != "all" and (feat.feature_type or "").lower() != feature_type:
            continue
        stats = extract_feature_stats(feat)
        if not stats or stats.get("log2fc") is None or stats.get("pval") is None:
            continue
        rows.append({"feature": feat.name, "log2FC": stats["log2fc"], "pvalue": stats["pval"]})

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["-log10p"] = -np.log10(df["pvalue"])
    return df


def render() -> None:
    st.header("Volcano Plot")
    st.caption("Differential expression: log2 fold-change vs -log10(p-value)")

    with db_session() as db:
        datasets = list_datasets_for_dropdown(db)

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

        fig = create_volcano_plot(
            df,
            fc_col="log2FC",
            pval_col="-log10p",
            fc_thresh=fc_thresh,
            p_thresh=p_thresh,
        )

        st.plotly_chart(fig, width='stretch')
        st.download_button(
            "Download CSV",
            data=df.to_csv(index=False),
            file_name="volcano_data.csv",
            mime="text/csv",
        )

