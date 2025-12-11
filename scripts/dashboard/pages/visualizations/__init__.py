from __future__ import annotations

import streamlit as st

from scripts.dashboard.pages.visualizations.convergence import render as render_convergence
from scripts.dashboard.pages.visualizations.dendrogram import render as render_dendrogram
from scripts.dashboard.pages.visualizations.heatmap import render as render_heatmap
from scripts.dashboard.pages.visualizations.network import render as render_network
from scripts.dashboard.pages.visualizations.pathway_bars import render as render_pathway_bars
from scripts.dashboard.pages.visualizations.pca import render as render_pca
from scripts.dashboard.pages.visualizations.program_signature import render as render_program_signature
from scripts.dashboard.pages.visualizations.volcano import render as render_volcano


def render_visualizations_page() -> None:
    st.title("📊 Interactive Visualizations")
    st.markdown("Explore datasets and signatures with interactive plots.")

    tabs = st.tabs([
        "Volcano",
        "Heatmap",
        "PCA",
        "Network",
        "Program-Signature",
        "Dataset Similarity",
        "Pathway Enrichment",
        "Cross-Omics Convergence",
    ])

    with tabs[0]:
        render_volcano()
    with tabs[1]:
        render_heatmap()
    with tabs[2]:
        render_pca()
    with tabs[3]:
        render_network()
    with tabs[4]:
        render_program_signature()
    with tabs[5]:
        render_dendrogram()
    with tabs[6]:
        render_pathway_bars()
    with tabs[7]:
        render_convergence()


__all__ = ["render_visualizations_page"]
