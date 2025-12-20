from __future__ import annotations

import math
from typing import List, Set

import pandas as pd
import plotly.express as px
import streamlit as st

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.database.models import Feature
from scripts.dashboard.db_session import db_session


def _run_enrichment(dataset_id: str, top_n: int = 20) -> pd.DataFrame:
    # Placeholder: re-use existing enrichment helper if available
    # Here we fetch features for the dataset and run pathway_enrichment
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(Feature.datasets)
            .filter_by(id=dataset_id)
            .all()
        )
    genes = [f.name for f in feats if f.name is not None]
    if not genes:
        return pd.DataFrame()
    input_features: Set[str] = {str(g) for g in genes if g}
    if not input_features:
        return pd.DataFrame()
    results = perform_pathway_enrichment(
        input_features=input_features,
        input_feature_types={"gene"},
    )
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    if "p_value" in df.columns:
        df["neg_log10_p"] = df["p_value"].apply(lambda p: -math.log10(p) if p > 0 else 0)
    return df.sort_values("p_value").head(top_n)


def render() -> None:
    st.header("Pathway Enrichment")
    st.caption("Horizontal bar chart of enriched pathways (−log10 p-value).")

    # Simple input: dataset ID text box
    dataset_id = st.text_input("Dataset ID (UUID)")
    top_n = st.slider("Top pathways", 5, 50, 20)

    if not dataset_id:
        st.info("Enter a Dataset ID to run enrichment.")
        return

    if st.button("Run Enrichment"):
        with st.spinner("Running pathway enrichment..."):
            df = _run_enrichment(dataset_id, top_n=top_n)

        if df.empty:
            st.warning("No enrichment results found.")
            return

        color_col = "database" if "database" in df.columns else None
        fig = px.bar(
            df,
            x="neg_log10_p" if "neg_log10_p" in df.columns else df.columns[0],
            y="pathway" if "pathway" in df.columns else df.columns[0],
            color=color_col,
            orientation="h",
            title="Pathway Enrichment",
        )
        fig.update_layout(yaxis={'categoryorder': 'total ascending'})
        st.plotly_chart(fig, use_container_width=True)

        st.download_button(
            "Download Results (CSV)",
            data=df.to_csv(index=False),
            file_name="pathway_enrichment.csv",
            mime="text/csv",
        )


__all__ = ["render"]
