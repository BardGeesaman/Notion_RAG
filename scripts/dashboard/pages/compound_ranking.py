"""Compound ranking dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 60) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 600) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_compound_ranking_page() -> None:
    """Render the compound ranking page."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Compound Ranking")
    st.caption("Multi-objective compound optimization and ranking.")

    tab1, tab2, tab3 = st.tabs(["Setup", "Ranked Table", "Pareto Plot"])

    with tab1:
        st.subheader("Compound Selection")
        
        # Compound ID input
        compound_ids_text = st.text_area(
            "Compound IDs (one per line)",
            placeholder="Enter compound UUIDs, one per line...",
            height=120,
            help="Enter compound UUIDs to rank. You can also fetch compounds from database queries."
        )
        
        # Parse compound IDs
        compound_ids = []
        if compound_ids_text.strip():
            compound_ids = [line.strip() for line in compound_ids_text.strip().split('\n') if line.strip()]
        
        st.write(f"**Compounds to rank:** {len(compound_ids)}")
        
        if len(compound_ids) > 500:
            st.error("Maximum 500 compounds allowed. Please reduce the number of compounds.")
        
        # Weight preset selector
        try:
            presets = _api_get("/api/ranking/presets")
            preset_options = [p["name"] for p in presets]
            preset_descriptions = {p["name"]: p["description"] for p in presets}
        except Exception as e:
            st.error(f"Failed to load presets: {e}")
            presets = []
            preset_options = ["balanced"]
            preset_descriptions = {"balanced": "Default balanced weights"}
        
        selected_preset = st.selectbox(
            "Weight Preset",
            options=preset_options,
            index=0,
            help="Choose a predefined weight preset or use custom weights below."
        )
        
        if selected_preset in preset_descriptions:
            st.info(f"**{selected_preset}:** {preset_descriptions[selected_preset]}")
        
        # Custom weights (collapsed by default)
        with st.expander("Custom Weights (Override Preset)", expanded=False):
            st.write("Adjust individual objective weights. Higher values = more importance.")
            
            # Get default weights from selected preset
            default_weights = {"potency": 0.35, "herg": 0.20, "alerts": 0.20, "logs": 0.15, "logp": 0.10}
            if presets:
                for preset in presets:
                    if preset["name"] == selected_preset:
                        default_weights = preset["weights"]
                        break
            
            col1, col2 = st.columns(2)
            
            with col1:
                potency_weight = st.slider(
                    "Potency Weight", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=default_weights.get("potency", 0.35),
                    step=0.05,
                    help="Weight for potency optimization (IC50/EC50)"
                )
                
                herg_weight = st.slider(
                    "hERG Safety Weight", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=default_weights.get("herg", 0.20),
                    step=0.05,
                    help="Weight for hERG cardiotoxicity safety"
                )
                
                alerts_weight = st.slider(
                    "Structural Alerts Weight", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=default_weights.get("alerts", 0.20),
                    step=0.05,
                    help="Weight for structural alert penalties (PAINS, etc.)"
                )
            
            with col2:
                logs_weight = st.slider(
                    "LogS (Solubility) Weight", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=default_weights.get("logs", 0.15),
                    step=0.05,
                    help="Weight for aqueous solubility"
                )
                
                logp_weight = st.slider(
                    "LogP (Lipophilicity) Weight", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=default_weights.get("logp", 0.10),
                    step=0.05,
                    help="Weight for lipophilicity optimization"
                )
            
            # Show total weight
            total_weight = potency_weight + herg_weight + alerts_weight + logs_weight + logp_weight
            st.write(f"**Total Weight:** {total_weight:.2f}")
            if abs(total_weight - 1.0) > 0.05:
                st.warning("Weights should sum to approximately 1.0 for best results.")
            
            custom_weights = {
                "potency": potency_weight,
                "herg": herg_weight,
                "alerts": alerts_weight,
                "logs": logs_weight,
                "logp": logp_weight
            }
            
            use_custom_weights = st.checkbox("Use Custom Weights", value=False)
        
        if not use_custom_weights:
            custom_weights = {}
        
        # Run ranking button
        col1, col2 = st.columns([1, 3])
        with col1:
            if st.button("Run Ranking", type="primary", disabled=len(compound_ids) == 0):
                if len(compound_ids) > 500:
                    st.error("Too many compounds. Maximum 500 allowed.")
                else:
                    with st.spinner("Running compound ranking..."):
                        try:
                            payload = {
                                "compound_ids": compound_ids,
                                "include_pareto": True
                            }
                            
                            if use_custom_weights:
                                payload["weights"] = custom_weights
                            else:
                                payload["preset"] = selected_preset
                            
                            results = _api_post("/api/ranking/score", payload, timeout=300)
                            st.session_state["ranking_results"] = results
                            st.success(f"Ranked {results['total_compounds']} compounds successfully!")
                            
                            if results["skipped_compounds"] > 0:
                                st.warning(f"Skipped {results['skipped_compounds']} compounds (missing potency data)")
                                
                        except Exception as e:
                            st.error(f"Ranking failed: {e}")
        
        with col2:
            if compound_ids:
                st.write(f"Ready to rank {len(compound_ids)} compounds")

    with tab2:
        st.subheader("Ranked Compounds")
        
        if "ranking_results" not in st.session_state:
            st.info("Run ranking in the Setup tab to see results here.")
        else:
            results = st.session_state["ranking_results"]
            rankings = results.get("rankings", [])
            
            if not rankings:
                st.warning("No compounds were successfully ranked.")
            else:
                # Build DataFrame
                data = []
                for ranking in rankings:
                    row = {
                        "Rank": ranking["rank"],
                        "Compound ID": ranking["compound_id"],
                        "SMILES": ranking["smiles"][:50] + "..." if len(ranking["smiles"]) > 50 else ranking["smiles"],
                        "Weighted Score": round(ranking["weighted_score"], 3),
                        "Pareto Rank": ranking["pareto_rank"]
                    }
                    
                    # Add objective scores
                    for obj in ranking["objectives"]:
                        obj_name = obj["name"].replace("_", " ").title()
                        row[f"{obj_name} (Raw)"] = obj["raw_value"]
                        row[f"{obj_name} (Norm)"] = round(obj["normalized"], 3)
                    
                    data.append(row)
                
                df = pd.DataFrame(data)
                
                # Display summary stats
                pareto_front = results.get("pareto_front", [])
                st.write(f"**Total Compounds:** {len(rankings)} | **Pareto Front:** {len(pareto_front)} compounds")
                
                # Display table with styling
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Weighted Score": st.column_config.ProgressColumn(
                            "Weighted Score",
                            help="Overall weighted score (0-1)",
                            min_value=0,
                            max_value=1,
                        ),
                        "Pareto Rank": st.column_config.NumberColumn(
                            "Pareto Rank",
                            help="Pareto dominance rank (1 = front)",
                            min_value=1,
                        ),
                    }
                )
                
                # Download button
                csv = df.to_csv(index=False)
                st.download_button(
                    label="Download Results as CSV",
                    data=csv,
                    file_name="compound_ranking_results.csv",
                    mime="text/csv"
                )

    with tab3:
        st.subheader("Pareto Plot")
        
        if "ranking_results" not in st.session_state:
            st.info("Run ranking in the Setup tab to see Pareto plot here.")
        else:
            results = st.session_state["ranking_results"]
            rankings = results.get("rankings", [])
            
            if not rankings:
                st.warning("No compounds were successfully ranked.")
            else:
                # Axis selectors
                col1, col2 = st.columns(2)
                
                # Get available objectives
                objectives = set()
                for ranking in rankings:
                    for obj in ranking["objectives"]:
                        objectives.add(obj["name"])
                
                objective_options = sorted(list(objectives))
                
                with col1:
                    x_objective = st.selectbox(
                        "X-Axis Objective",
                        options=objective_options,
                        index=objective_options.index("potency") if "potency" in objective_options else 0
                    )
                
                with col2:
                    y_options = objective_options + ["liability_aggregate"]
                    y_objective = st.selectbox(
                        "Y-Axis Objective",
                        options=y_options,
                        index=y_options.index("liability_aggregate") if "liability_aggregate" in y_options else 0
                    )
                
                # Build plot data
                plot_data = []
                for ranking in rankings:
                    obj_dict = {obj["name"]: obj["normalized"] for obj in ranking["objectives"]}
                    
                    x_val = obj_dict.get(x_objective, 0.0)
                    
                    if y_objective == "liability_aggregate":
                        # Aggregate liability (inverted - higher = worse)
                        herg_liability = 1.0 - obj_dict.get("herg", 0.0)
                        alerts_liability = 1.0 - obj_dict.get("alerts", 0.0)
                        logs_liability = 1.0 - obj_dict.get("logs", 0.0)
                        y_val = (herg_liability + alerts_liability + logs_liability) / 3.0
                    else:
                        y_val = obj_dict.get(y_objective, 0.0)
                    
                    plot_data.append({
                        "compound_id": ranking["compound_id"],
                        "smiles": ranking["smiles"],
                        "x_value": x_val,
                        "y_value": y_val,
                        "pareto_rank": ranking["pareto_rank"],
                        "weighted_score": ranking["weighted_score"],
                        "is_frontier": ranking["pareto_rank"] == 1
                    })
                
                plot_df = pd.DataFrame(plot_data)
                
                # Create scatter plot
                fig = px.scatter(
                    plot_df,
                    x="x_value",
                    y="y_value",
                    color="pareto_rank",
                    size="weighted_score",
                    hover_data=["compound_id", "smiles"],
                    title=f"Pareto Plot: {x_objective.title()} vs {y_objective.replace('_', ' ').title()}",
                    labels={
                        "x_value": x_objective.replace("_", " ").title(),
                        "y_value": y_objective.replace("_", " ").title(),
                        "pareto_rank": "Pareto Rank"
                    },
                    color_continuous_scale="Viridis_r"
                )
                
                # Highlight Pareto front
                frontier_df = plot_df[plot_df["is_frontier"]]
                if not frontier_df.empty:
                    # Sort frontier points by x_value for line drawing
                    frontier_sorted = frontier_df.sort_values("x_value")
                    
                    fig.add_trace(go.Scatter(
                        x=frontier_sorted["x_value"],
                        y=frontier_sorted["y_value"],
                        mode="lines+markers",
                        name="Pareto Front",
                        line=dict(color="red", width=2),
                        marker=dict(color="red", size=8, symbol="diamond")
                    ))
                
                fig.update_layout(
                    height=600,
                    showlegend=True,
                    hovermode="closest"
                )
                
                st.plotly_chart(fig, use_container_width=True)
                
                # Summary stats
                st.write(f"**Pareto Front:** {len(frontier_df)} compounds out of {len(plot_df)} total")
                
                if not frontier_df.empty:
                    st.write("**Pareto Front Compounds:**")
                    frontier_summary = frontier_df[["compound_id", "x_value", "y_value", "weighted_score"]].copy()
                    frontier_summary.columns = ["Compound ID", f"{x_objective.title()}", f"{y_objective.replace('_', ' ').title()}", "Weighted Score"]
                    st.dataframe(frontier_summary, hide_index=True, use_container_width=True)


__all__ = ["render_compound_ranking_page"]
