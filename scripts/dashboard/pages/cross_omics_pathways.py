"""Cross-Omics Pathway Analysis dashboard page."""
from __future__ import annotations

import pandas as pd
import plotly.express as px
import streamlit as st

from amprenta_rag.analysis.cross_omics_pathways import get_cross_omics_enrichment
from amprenta_rag.database.models import Program
from amprenta_rag.database.session import db_session


def render_cross_omics_pathways_page() -> None:
    st.header("🧬 Cross-Omics Pathway Analysis")
    st.caption("Identify pathways converging across transcriptomics, proteomics, metabolomics, and lipidomics.")

    with db_session() as db:
        programs = db.query(Program).order_by(Program.created_at.desc()).limit(100).all()

    if not programs:
        st.info("No programs found.")
        return

    prog_map = {p.name: p for p in programs}
    selected = st.selectbox("Select program", list(prog_map.keys()))
    program = prog_map[selected]

    # Omics toggles (currently informational; analysis uses all available)
    st.multiselect(
        "Omics to include",
        ["transcriptomics", "proteomics", "metabolomics", "lipidomics"],
        default=["transcriptomics", "proteomics", "metabolomics", "lipidomics"],
        disabled=True,
        help="All available omics are included automatically.",
    )

    try:
        result = get_cross_omics_enrichment(program.id)
    except Exception as exc:  # pragma: no cover - UI fallback
        st.error(f"Failed to run pathway analysis: {exc}")
        return

    pathways = result.pathways
    if not pathways:
        st.info("No convergent pathways found.")
        return

    # Convergence heatmap
    rows = []
    for p in pathways:
        for omics, feats in p.matched_by_omics.items():
            rows.append(
                {"Pathway": p.name, "Omics": omics, "Count": len(feats or [])}
            )
    df_heat = pd.DataFrame(rows)
    fig = px.imshow(
        df_heat.pivot_table(index="Pathway", columns="Omics", values="Count", fill_value=0),
        labels=dict(x="Omics", y="Pathway", color="Matched features"),
        aspect="auto",
        origin="lower",
        color_continuous_scale="Blues",
    )
    st.subheader("Convergence Heatmap")
    st.plotly_chart(fig, use_container_width=True)

    # Pathway table
    table_rows = []
    for p in pathways:
        table_rows.append(
            {
                "Pathway": p.name,
                "Source": p.source,
                "Adj P-value": p.adjusted_p_value,
                "Enrichment": p.enrichment_ratio,
                "Transcriptomics": len(p.matched_by_omics.get("transcriptomics", [])),
                "Proteomics": len(p.matched_by_omics.get("proteomics", [])),
                "Metabolomics": len(p.matched_by_omics.get("metabolomics", [])),
                "Lipidomics": len(p.matched_by_omics.get("lipidomics", [])),
            }
        )
    df_table = pd.DataFrame(table_rows)
    st.subheader("Pathways")
    st.dataframe(df_table, use_container_width=True, hide_index=True)

    # Feature list per pathway
    st.markdown("---")
    st.subheader("Pathway Feature Details")
    for p in pathways:
        with st.expander(f"{p.name} ({p.source})"):
            for omics in ["transcriptomics", "proteomics", "metabolomics", "lipidomics"]:
                feats = p.matched_by_omics.get(omics, [])
                st.markdown(f"**{omics.title()}** ({len(feats)}):")
                if feats:
                    st.write(", ".join(feats))
                else:
                    st.caption("None")

