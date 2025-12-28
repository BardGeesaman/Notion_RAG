"""HTS QC and Triage page."""
from __future__ import annotations

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from amprenta_rag.analysis.hts_qc import (
    detect_plate_format,
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

    # Campaign selector outside tabs
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

    # Create tabs
    tab1, tab2 = st.tabs(["QC Metrics", "Interactive Plate"])
    
    with tab1:
        _render_qc_metrics(campaign)
    
    with tab2:
        _render_interactive_plate(campaign)


def _render_qc_metrics(campaign) -> None:
    """Render the QC metrics tab with existing functionality."""
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

    # Bayesian Dose-Response Toggle
    st.markdown("---")
    
    use_bayesian = st.checkbox("Use Bayesian Dose-Response", value=False, 
                               help="Fit 4PL model with uncertainty quantification")

    if use_bayesian:
        with st.expander("Prior Configuration", expanded=False):
            col_a, col_b = st.columns(2)
            with col_a:
                ec50_mean = st.slider("EC50 Prior Mean (log)", -5.0, 5.0, 0.0, 0.1)
                ec50_sd = st.slider("EC50 Prior SD", 0.1, 3.0, 1.0, 0.1)
            with col_b:
                hill_mean = st.slider("Hill Prior Mean", -2.0, 4.0, 1.0, 0.1)
                hill_sd = st.slider("Hill Prior SD", 0.1, 3.0, 2.0, 0.1)
        
        # Store in session state for later use
        st.session_state["bayesian_prior_config"] = {
            "ec50_prior_mean": ec50_mean,
            "ec50_prior_sd": ec50_sd,
            "hill_prior_mean": hill_mean,
            "hill_prior_sd": hill_sd,
        }

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


def _render_interactive_plate(campaign) -> None:
    """Render the interactive plate visualization tab."""
    st.subheader("Interactive Plate Visualization")
    
    # Color mode selector
    color_mode = st.selectbox(
        "Color mode",
        options=["Normalized Value", "Z-Score", "Hit Status", "Raw Value"],
        key="plate_color_mode"
    )
    
    # Filter controls
    col1, col2 = st.columns(2)
    with col1:
        show_all = st.checkbox("Show all wells", value=True)
    with col2:
        z_threshold = st.slider(
            "Z-Score threshold", 
            min_value=-5.0, 
            max_value=5.0, 
            value=2.0, 
            step=0.1
        )
    
    # Fetch plate data
    try:
        wells = get_plate_heatmap_data(campaign.id, include_all=show_all)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to load plate data: {exc}")
        return
    
    # Apply Z-score threshold filter
    if wells and z_threshold is not None:
        wells = [w for w in wells if w.z_score is None or abs(w.z_score) >= z_threshold]
    
    if not wells:
        st.info("No well data available for the selected filters.")
        return
    
    # Build plate grid
    grid, customdata, max_row, max_col = _build_plate_grid(wells, color_mode)
    
    # Display plate format
    plate_format = detect_plate_format(max_row, max_col)
    st.caption(f"Plate format: {plate_format}")
    
    # Create Plotly heatmap
    row_labels = [chr(ord("A") + i) for i in range(max_row + 1)]
    col_labels = [str(i + 1).zfill(2) for i in range(max_col + 1)]
    
    # Set color scale based on mode
    if color_mode == "Normalized Value":
        colorscale = "RdYlGn_r"
        colorbar_title = "Normalized Value"
    elif color_mode == "Z-Score":
        colorscale = "RdBu_r"
        colorbar_title = "Z-Score"
    elif color_mode == "Hit Status":
        colorscale = [[0, "lightgray"], [1, "red"]]
        colorbar_title = "Hit Status"
    else:  # Raw Value
        colorscale = "Viridis"
        colorbar_title = "Raw Value"
    
    fig = go.Figure(data=go.Heatmap(
        z=grid,
        x=col_labels,
        y=row_labels,
        colorscale=colorscale,
        colorbar=dict(title=colorbar_title),
        customdata=customdata,
        hovertemplate="Well: %{y}%{x}<br>Value: %{z:.3f}<br>Compound: %{customdata[0]}<extra></extra>",
        showscale=True
    ))
    
    fig.update_layout(
        title=f"Interactive Plate - {color_mode}",
        xaxis_title="Column",
        yaxis_title="Row",
        width=800,
        height=600
    )
    
    # Display with selection capability
    selected_data = st.plotly_chart(
        fig, 
        use_container_width=True,
        on_select="rerun",
        selection_mode=["box", "points"]
    )
    
    # Well detail panel
    if selected_data and hasattr(selected_data, "selection"):
        selection = selected_data.selection
        if selection and (selection.points or selection.box):
            _render_well_detail_panel(wells, selection, campaign)


def _build_plate_grid(wells, color_mode):
    """Build 2D grid and customdata for plate visualization."""
    # Find plate dimensions
    max_row, max_col = -1, -1
    well_data = {}
    
    for w in wells:
        row_idx, col_idx = _parse_well_position(w.well_position or "")
        if row_idx >= 0 and col_idx >= 0:
            max_row = max(max_row, row_idx)
            max_col = max(max_col, col_idx)
            well_data[(row_idx, col_idx)] = w
    
    if max_row < 0 or max_col < 0:
        return [], [], 0, 0
    
    # Build grid based on color mode
    grid = [[None for _ in range(max_col + 1)] for _ in range(max_row + 1)]
    customdata = [[[""] for _ in range(max_col + 1)] for _ in range(max_row + 1)]
    
    for row in range(max_row + 1):
        for col in range(max_col + 1):
            well = well_data.get((row, col))
            if well:
                # Set value based on color mode
                if color_mode == "Normalized Value":
                    grid[row][col] = well.normalized_value
                elif color_mode == "Z-Score":
                    grid[row][col] = well.z_score
                elif color_mode == "Hit Status":
                    grid[row][col] = 1 if well.hit_flag else 0
                else:  # Raw Value
                    grid[row][col] = well.normalized_value  # Placeholder - adjust if raw_value available
                
                # Set custom data for hover
                compound_id = str(well.compound_id) if well.compound_id else "N/A"
                customdata[row][col] = [compound_id]
            else:
                grid[row][col] = None
                customdata[row][col] = [""]
    
    return grid, customdata, max_row, max_col


def _render_well_detail_panel(wells, selection, campaign):
    """Render well detail panel for selected wells."""
    st.subheader("Selected Wells")
    
    if not selection:
        st.info("Select wells on the plate to see details.")
        return
    
    # Extract selected well positions from SelectionState
    selected_wells = []
    points = []
    
    # Handle points selection
    if hasattr(selection, 'points') and selection.points:
        points = selection.points
    
    # Handle box selection (if any)
    if hasattr(selection, 'box') and selection.box:
        # Box selection could be implemented later if needed
        st.info("Box selection not yet implemented. Please use point selection.")
        return
    
    if not points:
        st.info("No wells selected. Click on wells to see details.")
        return
    
    for point in points:
        if "x" in point and "y" in point:
            col_label = point["x"]
            row_label = point["y"]
            well_position = f"{row_label}{col_label}"
            
            # Find matching well data
            for well in wells:
                if well.well_position == well_position:
                    selected_wells.append(well)
                    break
    
    if not selected_wells:
        st.info("No well data found for selected positions.")
        return
    
    # Create detail table
    detail_rows = []
    for well in selected_wells:
        detail_rows.append({
            "Well": well.well_position or "",
            "Compound ID": str(well.compound_id) if well.compound_id else "N/A",
            "Normalized": well.normalized_value,
            "Z-Score": well.z_score,
            "Hit": "Yes" if well.hit_flag else "No"
        })
    
    df_details = pd.DataFrame(detail_rows)
    st.dataframe(df_details, use_container_width=True, hide_index=True)
    
    # CSV download for selected wells
    if len(selected_wells) > 0:
        csv = df_details.to_csv(index=False)
        st.download_button(
            "📥 Download Selected Wells (CSV)",
            data=csv,
            file_name=f"{campaign.campaign_id}_selected_wells.csv",
            mime="text/csv",
        )

