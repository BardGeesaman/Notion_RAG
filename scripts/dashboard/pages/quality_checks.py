from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.analysis.quality_metrics import compute_quality_score
from amprenta_rag.database.models import Dataset
from scripts.dashboard.db_session import db_session


def _load_dataset_qc(db: Session):
    datasets = db.query(Dataset).order_by(Dataset.updated_at.desc()).limit(200).all()
    rows = []
    for ds in datasets:
        qc = compute_quality_score(ds)
        rows.append(
            {
                "id": str(ds.id),
                "name": ds.name,
                "omics_type": ds.omics_type or "",
                "score": qc["score"],
                "status": qc["status"],
                "issues": "; ".join(qc["issues"]) if qc["issues"] else "",
                "feature_count": qc["metrics"]["feature_count"],
                "stats_coverage_pct": qc["metrics"]["stats_coverage_pct"],
            }
        )
    return pd.DataFrame(rows)


def render_quality_checks_page() -> None:
    st.title("✅ Quality Checks")
    st.caption("Dataset quality scores with traffic-light status.")

    status_filter = st.selectbox(
        "Filter by status",
        ["all", "high", "medium", "low"],
        index=0,
        key="qc_status_filter",
    )

    if st.button("Refresh", key="qc_refresh"):
        st.rerun()

    with st.spinner("Computing quality scores..."):
        with db_session() as db:
            df = _load_dataset_qc(db)

    if df.empty:
        st.info("No datasets available.")
        return

    if status_filter != "all":
        df = df[df["status"] == status_filter]

    st.dataframe(
        df.sort_values("score", ascending=False),
        width='stretch',
        hide_index=True,
    )

    st.download_button(
        "Export CSV",
        data=df.to_csv(index=False),
        file_name="dataset_quality.csv",
        mime="text/csv",
        key="qc_export",
    )

    st.markdown("**Detail**")
    selection = st.selectbox(
        "Select dataset",
        df["name"].tolist(),
        key="qc_detail_select",
    )
    detail = df[df["name"] == selection].iloc[0].to_dict()
    st.json(detail, expanded=False)

