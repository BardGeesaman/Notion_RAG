"""Multi-omics latent factor integration dashboard page."""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

try:
    from upsetplot import plot as upset_plot
except ImportError:
    upset_plot = None


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 600) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _safe_json(text: str) -> Any:
    if not text.strip():
        return None
    return json.loads(text)


def render_multi_omics_integration_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Multi-Omics Integration")
    st.caption("Create and run multi-omics latent factor experiments (MOFA-style).")

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Setup + Run", "Variance", "Factor Scatter", "Loadings", "Alluvial", "UpSet"])

    experiments: List[Dict[str, Any]] = []
    try:
        experiments = list(_api_get("/api/multi-omics/experiments") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load experiments: {e}")
        experiments = []

    exp_options = [e.get("id") for e in experiments if e.get("id")]
    default_exp = st.session_state.get("multi_omics_selected_experiment")
    if default_exp not in exp_options:
        default_exp = exp_options[0] if exp_options else None

    with tab1:
        st.subheader("Existing experiments")
        if experiments:
            st.dataframe(pd.DataFrame(experiments), use_container_width=True, hide_index=True)
        else:
            st.info("No experiments found.")

        st.divider()
        st.subheader("Create experiment")
        name = st.text_input("Name", placeholder="My Multi-Omics Experiment")
        description = st.text_input("Description (optional)", placeholder="Short description")
        n_factors = st.number_input("n_factors", min_value=1, max_value=200, value=10, step=1)
        convergence_mode = st.selectbox("convergence_mode", options=["fast", "medium", "slow"], index=0)

        dataset_ids_text = st.text_area(
            "dataset_ids (JSON)",
            value='[{"omics_type":"transcriptomics","dataset_id":"<uuid>"},{"omics_type":"proteomics","dataset_id":"<uuid>"}]',
            height=120,
        )
        sample_mapping_text = st.text_area(
            "sample_mapping (JSON, optional)",
            value='[{"sample_id":"S1","transcriptomics":"RNA_A","proteomics":"PROT_A"}]',
            height=120,
        )

        c1, c2 = st.columns(2)
        with c1:
            if st.button("Create", type="primary"):
                try:
                    dataset_ids = _safe_json(dataset_ids_text)
                    sample_mapping = _safe_json(sample_mapping_text)
                    payload = {
                        "name": name,
                        "description": description or None,
                        "dataset_ids": dataset_ids,
                        "sample_mapping": sample_mapping,
                        "n_factors": int(n_factors),
                        "convergence_mode": str(convergence_mode),
                    }
                    out = _api_post("/api/multi-omics/experiments", payload, timeout=60)
                    st.success(f"Created experiment: {out.get('id')}")
                    st.session_state["multi_omics_selected_experiment"] = out.get("id")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Create failed: {e}")
        with c2:
            if st.button("Run selected"):
                eid = st.session_state.get("multi_omics_selected_experiment") or default_exp
                if not eid:
                    st.error("No experiment selected.")
                else:
                    try:
                        out = _api_post(f"/api/multi-omics/experiments/{eid}/run", {}, timeout=1800)
                        st.success(f"Run complete. Factors created: {out.get('factors_created')}")
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Run failed: {e}")

    with tab2:
        st.subheader("Variance explained per factor")
        eid = st.selectbox(
            "Experiment (variance)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_variance_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                if not factors:
                    st.info("No factors found (run experiment first).")
                else:
                    df = pd.DataFrame(factors)
                    df["total_variance"] = df["variance_explained"].apply(
                        lambda v: float(sum((v or {}).values())) if isinstance(v, dict) else None
                    )
                    fig = px.bar(df, x="factor_index", y="total_variance", title="Total variance explained (sum across views)")
                    st.plotly_chart(fig, use_container_width=True)
                    st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load variance: {e}")
        else:
            st.info("No experiments available.")

    with tab3:
        st.subheader("Factor scatter (Factor1 vs Factor2)")
        eid = st.selectbox(
            "Experiment (scatter)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_scatter_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                idxs = [int(f.get("factor_index")) for f in factors if f.get("factor_index") is not None]
                if len(idxs) < 2:
                    st.info("Need at least 2 factors (run experiment first).")
                else:
                    f1 = st.selectbox("X factor", options=idxs, index=0, key="multi_omics_f1")
                    f2 = st.selectbox("Y factor", options=idxs, index=1, key="multi_omics_f2")
                    s1 = pd.DataFrame(_api_get(f"/api/multi-omics/experiments/{eid}/factors/{f1}/scores") or [])
                    s2 = pd.DataFrame(_api_get(f"/api/multi-omics/experiments/{eid}/factors/{f2}/scores") or [])
                    if s1.empty or s2.empty:
                        st.info("No scores found (run experiment first).")
                    else:
                        s1 = s1.rename(columns={"score": "x"}).set_index("sample_id")
                        s2 = s2.rename(columns={"score": "y"}).set_index("sample_id")
                        merged = s1.join(s2, how="inner").reset_index()
                        fig = px.scatter(merged, x="x", y="y", hover_name="sample_id", title=f"Factor {f1} vs Factor {f2}")
                        st.plotly_chart(fig, use_container_width=True)
                        st.dataframe(merged, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load scatter: {e}")
        else:
            st.info("No experiments available.")

    with tab4:
        st.subheader("Top loadings")
        eid = st.selectbox(
            "Experiment (loadings)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_loadings_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                idxs = [int(f.get("factor_index")) for f in factors if f.get("factor_index") is not None]
                if not idxs:
                    st.info("No factors found (run experiment first).")
                else:
                    fi = st.selectbox("Factor", options=idxs, index=0, key="multi_omics_load_factor")
                    limit = st.slider("Limit", min_value=10, max_value=200, value=50, step=10)
                    rows = list(
                        _api_get(f"/api/multi-omics/experiments/{eid}/factors/{fi}/loadings?limit={int(limit)}") or []
                    )
                    if not rows:
                        st.info("No loadings found.")
                    else:
                        df = pd.DataFrame(rows)
                        st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load loadings: {e}")
        else:
            st.info("No experiments available.")

    with tab5:
        st.subheader("Alluvial (Dataset Flow)")
        
        # Get available datasets
        datasets: List[Dict[str, Any]] = []
        try:
            datasets = list(_api_get("/api/multi-omics-viz/datasets") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load datasets: {e}")
            datasets = []
        
        if not datasets:
            st.info("No datasets available.")
        else:
            # Dataset selection
            dataset_options = [{"id": d["id"], "label": f"{d['name']} ({d['omics_type']})"} for d in datasets]
            selected_indices = st.multiselect(
                "Select datasets",
                options=range(len(dataset_options)),
                format_func=lambda i: dataset_options[i]["label"],
                default=list(range(min(3, len(dataset_options)))),  # Select first 3 by default
                key="alluvial_datasets"
            )
            
            # Top N features slider
            top_n = st.slider("Top N features", min_value=10, max_value=100, value=50, step=10, key="alluvial_top_n")
            
            if selected_indices:
                selected_dataset_ids = [dataset_options[i]["id"] for i in selected_indices]
                
                if st.button("Generate Alluvial", type="primary", key="alluvial_generate"):
                    try:
                        payload = {"dataset_ids": selected_dataset_ids}
                        result = _api_post("/api/multi-omics-viz/alluvial", payload, timeout=120)
                        
                        # Build Plotly Sankey
                        nodes = result.get("nodes", [])
                        links = result.get("links", [])
                        
                        if not nodes or not links:
                            st.info("No data to visualize.")
                        else:
                            # Limit to top N features for better visualization
                            feature_nodes = [n for n in nodes if not any(d["id"] in n["id"] for d in datasets)]
                            if len(feature_nodes) > top_n:
                                # Keep top N features by total link value
                                feature_values = {}
                                for link in links:
                                    target_id = nodes[link["target"]]["id"]
                                    if target_id not in [d["id"] for d in datasets]:  # Is a feature
                                        feature_values[target_id] = feature_values.get(target_id, 0) + link["value"]
                                
                                top_features = sorted(feature_values.items(), key=lambda x: x[1], reverse=True)[:top_n]
                                top_feature_ids = {f[0] for f in top_features}
                                
                                # Filter nodes and links
                                filtered_nodes = [n for n in nodes if n["id"] in top_feature_ids or any(d["id"] in n["id"] for d in datasets)]
                                node_id_map = {old_id: new_id for new_id, old_id in enumerate([n["id"] for n in filtered_nodes])}
                                
                                filtered_links = []
                                for link in links:
                                    source_id = nodes[link["source"]]["id"]
                                    target_id = nodes[link["target"]]["id"]
                                    if source_id in node_id_map and target_id in node_id_map:
                                        filtered_links.append({
                                            "source": node_id_map[source_id],
                                            "target": node_id_map[target_id],
                                            "value": link["value"]
                                        })
                                
                                nodes = filtered_nodes
                                links = filtered_links
                                
                                # Update link indices
                                for i, link in enumerate(links):
                                    links[i] = {
                                        "source": link["source"],
                                        "target": link["target"], 
                                        "value": link["value"]
                                    }
                            
                            fig = go.Figure(data=[go.Sankey(
                                node=dict(
                                    pad=15,
                                    thickness=20,
                                    line=dict(color="black", width=0.5),
                                    label=[n["label"] for n in nodes],
                                    color=[n["color"] for n in nodes]
                                ),
                                link=dict(
                                    source=[link["source"] for link in links],
                                    target=[link["target"] for link in links],
                                    value=[link["value"] for link in links]
                                )
                            )])
                            
                            fig.update_layout(
                                title_text=f"Multi-Omics Data Flow (Top {len([n for n in nodes if not any(d['id'] in n['id'] for d in datasets)])} Features)",
                                font_size=10,
                                height=600
                            )
                            st.plotly_chart(fig, use_container_width=True)
                            
                            # Show summary stats
                            st.write(f"**Nodes:** {len(nodes)} | **Links:** {len(links)}")
                            
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Failed to generate alluvial: {e}")

    with tab6:
        st.subheader("UpSet Plot (Feature Overlaps)")
        
        # Reuse dataset selection
        if not datasets:
            st.info("No datasets available.")
        else:
            selected_indices = st.multiselect(
                "Select datasets",
                options=range(len(dataset_options)),
                format_func=lambda i: dataset_options[i]["label"],
                default=list(range(min(4, len(dataset_options)))),  # Select first 4 by default
                key="upset_datasets"
            )
            
            if selected_indices:
                selected_dataset_ids = [dataset_options[i]["id"] for i in selected_indices]
                
                if st.button("Generate UpSet Plot", type="primary", key="upset_generate"):
                    try:
                        payload = {"dataset_ids": selected_dataset_ids}
                        result = _api_post("/api/multi-omics-viz/upset", payload, timeout=120)
                        
                        sets_data = result.get("sets", [])
                        intersections = result.get("intersections", [])
                        matrix_items = result.get("matrix", [])
                        
                        if not sets_data or not intersections:
                            st.info("No overlap data to visualize.")
                        else:
                            if upset_plot is None:
                                st.warning("UpSet plotting requires: `pip install upsetplot`")
                                
                                # Fallback table view
                                st.subheader("Dataset Sizes")
                                sets_df = pd.DataFrame(sets_data)
                                st.dataframe(sets_df, use_container_width=True, hide_index=True)
                                
                                st.subheader("Top Intersections")
                                # Sort intersections by count
                                sorted_intersections = sorted(intersections, key=lambda x: x["count"], reverse=True)
                                intersect_df = pd.DataFrame(sorted_intersections[:20])  # Top 20
                                if not intersect_df.empty:
                                    intersect_df["dataset_names"] = intersect_df["sets"].apply(
                                        lambda set_ids: ", ".join([
                                            next((s["name"] for s in sets_data if s["id"] == sid), sid)
                                            for sid in set_ids
                                        ])
                                    )
                                st.dataframe(intersect_df, use_container_width=True, hide_index=True)
                                
                            else:
                                # Generate UpSet plot using upsetplot
                                import matplotlib.pyplot as plt
                                
                                # Convert to upsetplot format
                                set_names = [s["name"] for s in sets_data]
                                
                                # Build membership matrix
                                feature_memberships = {}
                                for item in matrix_items:
                                    feature = item["label"]
                                    presence = item["presence"]
                                    feature_memberships[feature] = {
                                        set_names[i]: bool(presence[i]) for i in range(len(set_names))
                                    }
                                
                                if feature_memberships:
                                    # Create DataFrame for upsetplot
                                    membership_df = pd.DataFrame.from_dict(feature_memberships, orient='index')
                                    
                                    # Generate plot
                                    fig = plt.figure(figsize=(12, 8))
                                    # Pass boolean DataFrame directly to upset_plot
                                    upset_plot(membership_df, 
                                              subset_size='count',
                                              show_counts=True,
                                              sort_by='cardinality',
                                              sort_categories_by='cardinality')
                                    plt.suptitle("Feature Overlaps Across Datasets", y=0.98)
                                    
                                    st.pyplot(fig, use_container_width=True)
                                    plt.close()
                                    
                                    # Show summary stats
                                    total_features = len(feature_memberships)
                                    total_intersections = len(intersections)
                                    st.write(f"**Total Features:** {total_features} | **Intersections:** {total_intersections}")
                                    
                                    # Show top intersections table
                                    st.subheader("Top Intersections")
                                    sorted_intersections = sorted(intersections, key=lambda x: x["count"], reverse=True)
                                    top_intersections = pd.DataFrame(sorted_intersections[:10])
                                    if not top_intersections.empty:
                                        top_intersections["dataset_names"] = top_intersections["sets"].apply(
                                            lambda set_ids: ", ".join([
                                                next((s["name"] for s in sets_data if s["id"] == sid), sid)
                                                for sid in set_ids
                                            ])
                                        )
                                    st.dataframe(top_intersections, use_container_width=True, hide_index=True)
                                
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Failed to generate UpSet plot: {e}")


__all__ = ["render_multi_omics_integration_page"]


