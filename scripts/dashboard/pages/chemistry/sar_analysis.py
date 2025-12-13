"""SAR analysis tab."""

from __future__ import annotations

import math

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


def render_sar_tab(
    tab,
    *,
    db_session,
    compound_model,
    get_compound_properties,
    calculate_lipinski,
    get_activity_data,
    detect_activity_cliffs,
    find_common_core,
    decompose_rgroups,
    get_rgroup_statistics,
) -> None:
    """Render the SAR Analysis tab."""
    with tab:
        st.subheader("SAR Analysis")
        props_df = get_compound_properties()
        if props_df.empty:
            st.info("No compounds available for SAR analysis.")
            return

        lipinski_results = props_df.apply(lambda r: calculate_lipinski(r.get("smiles", "")), axis=1)
        props_df["lipinski_compliant"] = [res.get("passes_ro5", False) for res in lipinski_results]
        compliant = int(props_df["lipinski_compliant"].sum())
        st.metric("Compounds", len(props_df))
        st.metric("Lipinski Compliant", compliant)

        col_a, col_b = st.columns(2)
        with col_a:
            fig_mw = px.histogram(props_df, x="molecular_weight", nbins=40, title="Molecular Weight")
            st.plotly_chart(fig_mw, use_container_width=True)
        with col_b:
            fig_logp = px.histogram(props_df, x="logp", nbins=40, title="LogP")
            st.plotly_chart(fig_logp, use_container_width=True)

        act_df = get_activity_data()
        if not act_df.empty:
            merged = act_df.merge(
                props_df[["compound_id", "logp"]],
                on="compound_id",
                how="left",
            )
            merged_clean = merged.dropna(subset=["molecular_weight", "logp"])
            if not merged_clean.empty:
                scatter = px.scatter(
                    merged_clean,
                    x="logp",
                    y="molecular_weight",
                    title="LogP vs Molecular Weight (HTS)",
                    hover_data=["compound_id", "smiles"],
                )
                st.plotly_chart(scatter, use_container_width=True)

            cliff_thresh = st.slider("Cliff similarity threshold", min_value=0.5, max_value=1.0, value=0.8, step=0.01)
            cliffs = detect_activity_cliffs(similarity_threshold=cliff_thresh)
            if cliffs:
                st.markdown("**Activity Cliffs (top 50 by activity difference)**")
                st.dataframe(cliffs[:50], hide_index=True, use_container_width=True)
            else:
                st.info("No activity cliffs detected at the selected threshold.")
        else:
            st.info("No HTS activity data available for cliff detection.")

        st.markdown("### Activity Cliffs Network")
        cliffs = detect_activity_cliffs()
        if cliffs:
            nodes = {}
            for c in cliffs:
                nodes[c["compound_1"]] = nodes.get(c["compound_1"], {"id": c["compound_1"]})
                nodes[c["compound_2"]] = nodes.get(c["compound_2"], {"id": c["compound_2"]})

            ids = list(nodes.keys())
            n = len(ids)
            for idx, cid in enumerate(ids):
                angle = 2 * math.pi * idx / max(n, 1)
                nodes[cid]["x"] = math.cos(angle)
                nodes[cid]["y"] = math.sin(angle)

            edge_x, edge_y, edge_text = [], [], []
            for c in cliffs:
                a, b = c["compound_1"], c["compound_2"]
                x0, y0 = nodes[a]["x"], nodes[a]["y"]
                x1, y1 = nodes[b]["x"], nodes[b]["y"]
                edge_x += [x0, x1, None]
                edge_y += [y0, y1, None]
                edge_text.append(
                    f"{a} ↔ {b}<br>similarity={c['similarity']:.2f}<br>fold_change={c['fold_change']:.1f}"
                )

            edge_trace = go.Scatter(
                x=edge_x,
                y=edge_y,
                line=dict(width=1, color="lightgray"),
                hoverinfo="text",
                text=edge_text,
                mode="lines",
            )

            node_x = [nodes[cid]["x"] for cid in ids]
            node_y = [nodes[cid]["y"] for cid in ids]
            node_text = ids

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                text=node_text,
                textposition="bottom center",
                marker=dict(size=12, color="tomato"),
                hoverinfo="text",
            )

            fig = go.Figure(data=[edge_trace, node_trace])
            fig.update_layout(
                showlegend=False,
                xaxis=dict(showgrid=False, zeroline=False, visible=False),
                yaxis=dict(showgrid=False, zeroline=False, visible=False),
                margin=dict(l=10, r=10, t=10, b=10),
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No activity cliffs detected.")

        st.markdown("### R-Group Decomposition")
        st.markdown("Find common core and decompose compounds into R-groups for SAR analysis.")

        with db_session() as db:
            compounds = db.query(compound_model).limit(100).all()
            if compounds:
                compound_options = {f"{c.compound_id} - {c.smiles[:30]}": c.smiles for c in compounds}
                selected_compounds = st.multiselect(
                    "Select compounds for R-group analysis",
                    list(compound_options.keys()),
                    help="Select 2 or more compounds to find common core",
                )

                if len(selected_compounds) >= 2:
                    selected_smiles = [compound_options[c] for c in selected_compounds]

                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button("Find Common Core", type="primary"):
                            with st.spinner("Finding maximum common substructure..."):
                                core_smarts = find_common_core(selected_smiles)
                                if core_smarts:
                                    st.session_state["rgroup_core"] = core_smarts
                                    st.success(f"Common core found: {core_smarts}")
                                else:
                                    st.error("Could not find common core. Try different compounds.")

                    if st.session_state.get("rgroup_core"):
                        core_smarts = st.session_state["rgroup_core"]
                        st.markdown(f"**Core SMARTS:** `{core_smarts}`")

                        with col2:
                            if st.button("Decompose R-Groups", type="primary"):
                                with st.spinner("Decomposing compounds..."):
                                    decomposition = decompose_rgroups(selected_smiles, core_smarts)
                                    if decomposition:
                                        st.session_state["rgroup_decomp"] = decomposition
                                        st.success(f"Decomposed {len(decomposition)} compounds")
                                    else:
                                        st.error("R-group decomposition failed")

                        if st.session_state.get("rgroup_decomp"):
                            decomp = st.session_state["rgroup_decomp"]

                            st.markdown("**Decomposition Results:**")
                            decomp_df = pd.DataFrame(decomp)
                            st.dataframe(decomp_df, hide_index=True, use_container_width=True)

                            stats = get_rgroup_statistics(decomp)
                            if stats:
                                st.markdown("**R-Group Statistics:**")
                                for position, counts in stats.items():
                                    st.markdown(f"**{position}:**")
                                    stats_df = pd.DataFrame(
                                        [
                                            {"R-Group": smiles, "Frequency": count}
                                            for smiles, count in sorted(counts.items(), key=lambda x: x[1], reverse=True)
                                        ]
                                    )
                                    st.dataframe(stats_df, hide_index=True, use_container_width=True)
                else:
                    st.info("Select at least 2 compounds to perform R-group analysis.")
            else:
                st.info("No compounds available for R-group analysis.")

