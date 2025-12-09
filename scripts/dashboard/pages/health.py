"""Health Check page for the Streamlit dashboard."""

from __future__ import annotations

from datetime import datetime

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from amprenta_rag.signatures.signature_validation import validate_all_signatures
from scripts.dashboard.db_session import db_session


def infer_status(ds):
    # Use ingestion_status (actual field) instead of qc_status (doesn't exist)
    ingestion_status = getattr(ds, "ingestion_status", None)
    
    if ingestion_status == "complete" and ds.features and len(ds.features) > 0:
        return "complete"
    if ingestion_status == "failed" or not ds.features:
        return "incomplete"
    if ingestion_status == "pending" or (hasattr(ds, "created_at") and (datetime.utcnow() - ds.created_at).total_seconds() < 1800):
        return "pending"
    if hasattr(ds, "updated_at") and (datetime.utcnow() - ds.updated_at).total_seconds() > 6 * 3600:
        return "stuck"
    return "in_progress"


def render_health_page() -> None:
    """
    Render the Health Check page showing system status and metrics.

    Features:
    - Database connection status
    - System health metrics
    - Data quality indicators
    """
    st.header("🏥 System Health Dashboard")
    st.info(
        "This page shows ingestion/back-pressure status, QC summaries, and signature validation at a glance. Statuses are inferred using the same rules as `scripts/list_ingestion_status.py`."
    )

    # All datasets and calculated statuses
    with db_session() as db:
        datasets = db.query(Dataset).all()
        now = datetime.utcnow()
        rows = []
        status_counts = {k: 0 for k in ["complete", "pending", "in_progress", "incomplete", "stuck"]}
        for ds in datasets:
            s = infer_status(ds)
            status_counts[s] += 1
            age_h = (now - ds.created_at).total_seconds() / 3600 if hasattr(ds, "created_at") else 0
            upd_h = (now - ds.updated_at).total_seconds() / 3600 if hasattr(ds, "updated_at") else 0
            rows.append(
                {
                    "ID": str(ds.id),
                    "Name": ds.name,
                    "Omics": ds.omics_type,
                    "Status": s,
                    "Ingestion": getattr(ds, "ingestion_status", "unknown"),
                    "Age (h)": int(age_h),
                    "Updated (h ago)": int(upd_h),
                    "Source": getattr(ds, "external_ids", getattr(ds, "file_paths", "")),
                    "Num Features": len(ds.features) if ds.features else 0,
                }
            )
        total = len(datasets)
        pct_comp = 100 * status_counts["complete"] / total if total else 0
        100 * status_counts["pending"] / total if total else 0
        stuck = [r for r in rows if r["Status"] == "stuck"]
        # Summary metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Datasets", total)
        with col2:
            st.metric("% Complete", f"{pct_comp:.1f}%")
        with col3:
            st.metric("# Stuck > 6h", len(stuck))
    # Bar of status breakdown
    st.bar_chart(pd.DataFrame([status_counts]))

    # Table for stuck datasets and filterable drilldown
    st.subheader("Stuck/Filtered Datasets")
    filter_status = st.selectbox("Status", ["All"] + list(status_counts.keys()))
    filter_omics = st.selectbox("Omics", ["All"] + sorted(set(r["Omics"] for r in rows)))
    filter_qc = st.selectbox("QC", ["All"] + sorted(set(r["QC"] for r in rows if r["QC"])))
    table = rows
    if filter_status != "All":
        table = [r for r in table if r["Status"] == filter_status]
    if filter_omics != "All":
        table = [r for r in table if r["Omics"] == filter_omics]
    if filter_qc != "All":
        table = [r for r in table if r["QC"] == filter_qc]
    df = pd.DataFrame(table)
    if len(df) > 0:
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.info("No datasets to display with this filter.")

    # QC distribution
    st.subheader("QC Status Distribution")
    qc_dist = pd.Series([r["QC"] for r in rows if r["QC"]]).value_counts()
    st.bar_chart(qc_dist)

    # Table of worst signatures by validation
    st.subheader("Signatures With Lowest Evidence")
    all_vals = validate_all_signatures()
    all_vals = sorted(all_vals, key=lambda v: v.metrics.coverage or 0)[:5]
    if all_vals:
        t = {
            "Signature": [str(v.signature_id) for v in all_vals],
            "Coverage": [v.metrics.coverage for v in all_vals],
            "Matched": [v.metrics.num_matched_datasets for v in all_vals],
            "Total": [v.metrics.num_total_datasets for v in all_vals],
            "Mean Score": [v.metrics.mean_score for v in all_vals],
        }
        st.dataframe(pd.DataFrame(t))
    else:
        st.info("No weak/low-evidence signatures found.")

    # Instructions/help
    st.markdown("---")
    with st.expander("How to Read this Page"):
        st.markdown(
            """
            - **Statuses**: `complete` = QC PASS + has features; `pending` = ingested <30min ago; `stuck` = no updates >6h; others in progress/incomplete.
            - **QC**: Derived from missingness/features per dataset.
            - **Signature Validation**: Validations based on your empirical overlap/cross-dataset analysis.
            - **Use cases:** Spot backlog (stuck), areas needing data/QC action, and weak signatures in analytics.
            - See also the CLI tools: `list_ingestion_status.py`, `check_pipeline_health.py`.
            """
        )
