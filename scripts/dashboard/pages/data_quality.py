"""Data Quality Watcher page."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.analysis.quality_watcher import (
    get_low_quality_datasets,
    get_quality_summary,
    scan_all_datasets,
)


def render_data_quality_page() -> None:
    """Render the Data Quality Watcher page."""
    st.header("✅ Data Quality Watcher")
    st.caption("Monitor dataset quality scores, distribution, and low-quality alerts.")

    try:
        summary = get_quality_summary()
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load quality summary: {exc}")
        return

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Datasets", summary.get("total", 0))
    col2.metric("High", summary.get("high", 0))
    col3.metric("Medium", summary.get("medium", 0))
    col4.metric("Low", summary.get("low", 0), delta_color="inverse")

    st.markdown("---")

    dist_df = pd.DataFrame(
        {
            "Status": ["high", "medium", "low"],
            "Count": [
                summary.get("high", 0),
                summary.get("medium", 0),
                summary.get("low", 0),
            ],
        }
    )
    st.subheader("Quality Distribution")
    st.bar_chart(dist_df.set_index("Status"))

    st.markdown("---")
    st.subheader("Low-Quality Alerts")
    threshold = st.slider("Alert threshold", min_value=0, max_value=100, value=50, step=5)
    try:
        low_quality = get_low_quality_datasets(threshold=threshold)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load low-quality datasets: {exc}")
        return

    if low_quality:
        alerts = []
        for r in low_quality:
            alerts.append(
                {
                    "Dataset": r.dataset_name,
                    "Score": r.score,
                    "Status": r.status,
                    "Issues": "; ".join(r.issues) if r.issues else "",
                    "Feature Count": r.metrics.get("feature_count"),
                    "Stats Coverage (%)": r.metrics.get("stats_coverage_pct"),
                }
            )
        alerts_df = pd.DataFrame(alerts)
        st.dataframe(alerts_df, use_container_width=True, hide_index=True)
    else:
        st.success("No low-quality datasets under the current threshold.")

    st.markdown("---")
    st.subheader("All Datasets")
    try:
        reports = scan_all_datasets()
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load dataset quality reports: {exc}")
        return
    if reports:
        rows = []
        for r in reports:
            rows.append(
                {
                    "Dataset": r.dataset_name,
                    "Score": r.score,
                    "Status": r.status,
                    "Issues": "; ".join(r.issues) if r.issues else "",
                    "Feature Count": r.metrics.get("feature_count"),
                    "Stats Coverage (%)": r.metrics.get("stats_coverage_pct"),
                }
            )
        all_df = pd.DataFrame(rows)
        st.dataframe(all_df, use_container_width=True, hide_index=True)
    else:
        st.info("No datasets available to score.")
