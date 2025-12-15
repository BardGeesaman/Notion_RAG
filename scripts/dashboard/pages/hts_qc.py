"""HTS QC and Triage page."""
from __future__ import annotations

import pandas as pd
import plotly.express as px
import streamlit as st

from amprenta_rag.analysis.hts_qc import (
    get_hit_compounds,
    get_plate_heatmap_data,
    get_plate_qc_summary,
)
from amprenta_rag.database.models import HTSCampaign
from amprenta_rag.database.session import db_session


def _parse_well_position(pos: str) -> tuple[int, int]:
    """Convert well position like A01 to (row_index, col_index)."""
    if not pos:
        return (-1, -1)
    row = pos[0].upper()
    try:
        col = int(pos[1:])
    except Exception:
        col = -1
    row_idx = ord(row) - ord("A")
    col_idx = col - 1 if col > 0 else -1
    return row_idx, col_idx


def _make_heatmap_df(wells):
    rows = []
    for w in wells:
        row_idx, col_idx = _parse_well_position(w.well_position or "")
        rows.append(
            {
                "row_idx": row_idx,
                "col_idx": col_idx,
                "Well": w.well_position,
                "Value": w.normalized_value,
                "Hit": w.hit_flag,
            }
        )
    df = pd.DataFrame(rows)
    return df


def render_hts_qc_page() -> None:
    st.header("🧪 HTS QC & Triage")
    st.caption("Review plate QC, hit rates, and triage hits for HTS campaigns.")

    with db_session() as db:
        campaigns = (
            db.query(HTSCampaign).order_by(HTSCampaign.created_at.desc()).limit(100).all()
        )

    if not campaigns:
        st.info("No HTS campaigns found.")
        return

    camp_map = {f"{c.campaign_name} ({c.campaign_id})": c for c in campaigns}
    selection = st.selectbox("Select campaign", list(camp_map.keys()))
    campaign = camp_map[selection]

    try:
        summary = get_plate_qc_summary(campaign.id)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load QC summary: {exc}")
        return

    col1, col2, col3 = st.columns(3)
    col1.metric("Z' factor", f"{summary.z_prime:.3f}" if summary.z_prime is not None else "N/A")
    col2.metric("Hit rate", f"{summary.hit_rate:.2f}%")
    col3.metric("Hits", summary.hits)
    st.caption(f"Controls — Pos: {summary.pos_controls}, Neg: {summary.neg_controls}")

    st.markdown("---")
    st.subheader("Plate Heatmap")
    try:
        wells = get_plate_heatmap_data(campaign.id)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load plate data: {exc}")
        return

    if wells:
        df_heat = _make_heatmap_df(wells)
        if not df_heat.empty:
            max_row = df_heat["row_idx"].max()
            max_col = df_heat["col_idx"].max()
            # Build grid; fill missing with NaN
            grid = [[None for _ in range(max_col + 1)] for _ in range(max_row + 1)]
            for _, r in df_heat.iterrows():
                if r.row_idx >= 0 and r.col_idx >= 0:
                    grid[r.row_idx][r.col_idx] = r.Value
            row_labels = [chr(ord("A") + i) for i in range(max_row + 1)]
            col_labels = [str(i + 1).zfill(2) for i in range(max_col + 1)]
            fig = px.imshow(
                grid,
                labels=dict(x="Column", y="Row", color="Normalized Value"),
                x=col_labels,
                y=row_labels,
                aspect="auto",
                origin="upper",
            )
            st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("No well data available.")

    st.markdown("---")
    st.subheader("Hit Compounds")
    try:
        hits = get_hit_compounds(campaign.id)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load hits: {exc}")
        return

    if hits:
        rows = []
        for h in hits:
            rows.append(
                {
                    "Result ID": h.result_id,
                    "Compound ID": str(h.compound_id),
                    "Well": h.well_position or "",
                    "Normalized": h.normalized_value,
                    "Z-score": h.z_score,
                }
            )
        df_hits = pd.DataFrame(rows)
        st.dataframe(df_hits, use_container_width=True, hide_index=True)
        csv = df_hits.to_csv(index=False)
        st.download_button(
            "📥 Download Hits (CSV)",
            data=csv,
            file_name=f"{campaign.campaign_id}_hits.csv",
            mime="text/csv",
        )
    else:
        st.success("No hits flagged for this campaign.")

