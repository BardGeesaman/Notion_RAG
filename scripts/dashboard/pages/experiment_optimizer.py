"""Experiment Optimizer page (Bayesian optimization recommendations)."""

from __future__ import annotations

import json
from typing import Any, Dict, List

import streamlit as st
import plotly.graph_objects as go

from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds, recommend_multi_objective


DEFAULT_TESTED = [
    {"features": [0.0, 0.0], "activity": 0.1, "potency": 0.1, "herg": 0.9, "logs": -2.0},
    {"features": [0.5, 0.5], "activity": 0.6, "potency": 0.6, "herg": 0.5, "logs": -1.5},
    {"features": [1.0, 0.2], "activity": 0.8, "potency": 0.8, "herg": 0.3, "logs": -1.0},
    {"features": [0.2, 0.9], "activity": 0.4, "potency": 0.4, "herg": 0.7, "logs": -2.5},
    {"features": [0.8, 0.8], "activity": 0.9, "potency": 0.9, "herg": 0.2, "logs": -0.8},
]

DEFAULT_POOL = [
    {"id": "c1", "features": [0.5, 0.5]},
    {"id": "c2", "features": [1.5, 1.2]},
    {"id": "c3", "features": [2.2, 1.9]},
    {"id": "c4", "features": [3.0, 3.0]},
]


def _parse_json(text: str) -> List[Dict[str, Any]]:
    obj = json.loads(text)
    if not isinstance(obj, list):
        raise ValueError("Expected a JSON list")
    return obj


def _render_pareto_2d(tested, recommendations, objectives, pareto_front):
    """Render 2D Pareto plot for 2 objectives."""
    import pandas as pd
    
    # Extract objective values from tested compounds
    tested_x = [c.get(objectives[0], 0) for c in tested]
    tested_y = [c.get(objectives[1], 0) for c in tested]
    
    # Get Pareto front points
    pareto_indices = pareto_front.get("indices", [])
    pareto_x = [tested[i].get(objectives[0], 0) for i in pareto_indices if i < len(tested)]
    pareto_y = [tested[i].get(objectives[1], 0) for i in pareto_indices if i < len(tested)]
    
    # Sort Pareto points for line
    if pareto_x and pareto_y:
        pareto_points = sorted(zip(pareto_x, pareto_y))
        pareto_x_sorted, pareto_y_sorted = zip(*pareto_points)
    else:
        pareto_x_sorted, pareto_y_sorted = [], []
    
    # Extract recommendations predictions
    rec_x = [r["pred_means"].get(objectives[0], 0) for r in recommendations]
    rec_y = [r["pred_means"].get(objectives[1], 0) for r in recommendations]
    
    fig = go.Figure()
    
    # Tested compounds
    fig.add_trace(go.Scatter(
        x=tested_x, y=tested_y,
        mode='markers',
        marker=dict(size=10, color='blue', symbol='circle'),
        name='Tested',
        hovertemplate=f'{objectives[0]}: %{{x:.3f}}<br>{objectives[1]}: %{{y:.3f}}<extra></extra>'
    ))
    
    # Pareto frontier
    if pareto_x_sorted:
        fig.add_trace(go.Scatter(
            x=pareto_x_sorted, y=pareto_y_sorted,
            mode='lines',
            line=dict(color='red', dash='dash', width=2),
            name='Pareto Front',
            hoverinfo='skip'
        ))
    
    # Recommendations
    if rec_x:
        fig.add_trace(go.Scatter(
            x=rec_x, y=rec_y,
            mode='markers',
            marker=dict(size=14, color='red', symbol='star'),
            name='Recommended',
            hovertemplate=f'{objectives[0]}: %{{x:.3f}}<br>{objectives[1]}: %{{y:.3f}}<extra></extra>'
        ))
    
    fig.update_layout(
        title='Pareto Front Visualization',
        xaxis_title=objectives[0],
        yaxis_title=objectives[1],
        hovermode='closest',
        height=500
    )
    
    st.plotly_chart(fig, use_container_width=True)


def _render_pareto_3d(tested, recommendations, objectives, pareto_front):
    """Render 3D Pareto plot for 3 objectives."""
    # Extract objective values from tested compounds
    tested_x = [c.get(objectives[0], 0) for c in tested]
    tested_y = [c.get(objectives[1], 0) for c in tested]
    tested_z = [c.get(objectives[2], 0) for c in tested]
    
    # Mark Pareto front points
    pareto_indices = set(pareto_front.get("indices", []))
    tested_pareto = ['Pareto' if i in pareto_indices else 'Non-Pareto' for i in range(len(tested))]
    
    # Extract recommendations predictions
    rec_x = [r["pred_means"].get(objectives[0], 0) for r in recommendations]
    rec_y = [r["pred_means"].get(objectives[1], 0) for r in recommendations]
    rec_z = [r["pred_means"].get(objectives[2], 0) for r in recommendations]
    
    fig = go.Figure()
    
    # Tested compounds
    fig.add_trace(go.Scatter3d(
        x=tested_x, y=tested_y, z=tested_z,
        mode='markers',
        marker=dict(
            size=6,
            color=[1 if p == 'Pareto' else 0 for p in tested_pareto],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Pareto"),
        ),
        text=tested_pareto,
        name='Tested',
        hovertemplate=f'{objectives[0]}: %{{x:.3f}}<br>{objectives[1]}: %{{y:.3f}}<br>{objectives[2]}: %{{z:.3f}}<extra></extra>'
    ))
    
    # Recommendations
    if rec_x:
        fig.add_trace(go.Scatter3d(
            x=rec_x, y=rec_y, z=rec_z,
            mode='markers',
            marker=dict(size=8, color='red', symbol='diamond'),
            name='Recommended',
            hovertemplate=f'{objectives[0]}: %{{x:.3f}}<br>{objectives[1]}: %{{y:.3f}}<br>{objectives[2]}: %{{z:.3f}}<extra></extra>'
        ))
    
    fig.update_layout(
        title='3D Pareto Front Visualization',
        scene=dict(
            xaxis_title=objectives[0],
            yaxis_title=objectives[1],
            zaxis_title=objectives[2],
        ),
        height=600
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_experiment_optimizer_page() -> None:
    """Render the Experiment Optimizer page."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧠 Experiment Optimizer")
    st.caption("BoTorch-based Bayesian optimization to recommend next compounds to test.")

    # Optimization mode selection
    opt_mode = st.radio(
        "Optimization Mode",
        ["Single-Objective", "Multi-Objective"],
        horizontal=True,
    )
    
    # Multi-objective configuration
    objectives = ["potency"]
    objective_directions = ["maximize"]
    
    if opt_mode == "Multi-Objective":
        st.markdown("#### Multi-Objective Configuration")
        available_objectives = ["potency", "herg", "logs", "selectivity"]
        objectives = st.multiselect(
            "Select Objectives (minimum 2)",
            available_objectives,
            default=["potency", "herg"],
        )
        
        if len(objectives) < 2:
            st.warning("⚠️ Multi-objective mode requires at least 2 objectives")
        
        if objectives:
            st.markdown("**Objective Directions:**")
            objective_directions = []
            cols = st.columns(len(objectives))
            for i, obj in enumerate(objectives):
                with cols[i]:
                    direction = st.selectbox(
                        f"{obj}",
                        ["maximize", "minimize"],
                        index=0 if obj == "potency" else 1,
                        key=f"dir_{obj}"
                    )
                    objective_directions.append(direction)

    st.markdown("### Inputs")
    col1, col2 = st.columns(2)
    with col1:
        tested_json = st.text_area(
            "Tested compounds (JSON list of {features, activity})",
            value=json.dumps(DEFAULT_TESTED, indent=2),
            height=240,
        )
    with col2:
        pool_json = st.text_area(
            "Candidate pool (JSON list of {id, features})",
            value=json.dumps(DEFAULT_POOL, indent=2),
            height=240,
        )

    batch_size = st.number_input("Batch size", min_value=1, max_value=100, value=10, step=1)

    if st.button("Recommend", type="primary", use_container_width=True):
        try:
            tested = _parse_json(tested_json)
            pool = _parse_json(pool_json)
            
            if opt_mode == "Single-Objective":
                recs = recommend_next_compounds(
                    tested_compounds=tested,
                    candidate_pool=pool,
                    batch_size=int(batch_size),
                    objectives=["potency"],
                )
                st.session_state["optimizer_recs"] = recs
                st.session_state["optimizer_mode"] = "single"
                st.session_state["optimizer_tested"] = tested
            else:
                # Multi-objective mode
                if len(objectives) < 2:
                    st.error("Multi-objective mode requires at least 2 objectives")
                    return
                
                # Validate objectives exist in tested data
                missing_objs = []
                for obj in objectives:
                    if not all(obj in c for c in tested):
                        missing_objs.append(obj)
                
                if missing_objs:
                    st.error(f"Missing objectives in tested data: {', '.join(missing_objs)}")
                    st.info("Update your tested compounds JSON to include all selected objectives.")
                    return
                
                result = recommend_multi_objective(
                    tested_compounds=tested,
                    candidate_pool=pool,
                    objectives=objectives,
                    objective_directions=objective_directions,
                    batch_size=int(batch_size),
                )
                st.session_state["optimizer_recs"] = result["recommendations"]
                st.session_state["optimizer_pareto"] = result["pareto_front"]
                st.session_state["optimizer_mode"] = "multi"
                st.session_state["optimizer_objectives"] = objectives
                st.session_state["optimizer_tested"] = tested
                
        except ImportError as e:
            st.error(str(e))
            st.info("Install dependencies and restart: botorch, ax-platform, torch.")
            return
        except Exception as e:  # noqa: BLE001
            st.error(f"Recommendation failed: {e}")
            import traceback
            st.code(traceback.format_exc())
            return

    recs = st.session_state.get("optimizer_recs") or []
    mode = st.session_state.get("optimizer_mode", "single")
    
    st.markdown("### Recommendations")
    if not recs:
        st.info("Run Recommend to see ranked suggestions.")
        return

    st.dataframe(recs, use_container_width=True, hide_index=True)
    
    # Show Pareto visualization for multi-objective
    if mode == "multi":
        pareto_front = st.session_state.get("optimizer_pareto", {})
        objs = st.session_state.get("optimizer_objectives", [])
        tested = st.session_state.get("optimizer_tested", [])
        
        if pareto_front and objs and tested:
            st.markdown("### Pareto Front Visualization")
            n_pareto = len(pareto_front.get("indices", []))
            st.info(f"📊 {n_pareto} compounds on Pareto front")
            
            if len(objs) == 2:
                _render_pareto_2d(tested, recs, objs, pareto_front)
            elif len(objs) == 3:
                _render_pareto_3d(tested, recs, objs, pareto_front)
            else:
                st.info(f"Pareto visualization available for 2 or 3 objectives (you selected {len(objs)})")


__all__ = ["render_experiment_optimizer_page"]


