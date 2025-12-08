from __future__ import annotations

import streamlit as st

from scripts.dashboard.pages.visualizations.heatmap import render as render_heatmap
from scripts.dashboard.pages.visualizations.network import render as render_network
from scripts.dashboard.pages.visualizations.pca import render as render_pca
from scripts.dashboard.pages.visualizations.volcano import render as render_volcano


def render_visualizations_page() -> None:
    st.title("📊 Interactive Visualizations")
    st.markdown("Explore datasets and signatures with interactive plots.")

    tabs = st.tabs(["Volcano", "Heatmap", "PCA", "Network"])

    with tabs[0]:
        render_volcano()
    with tabs[1]:
        render_heatmap()
    with tabs[2]:
        render_pca()
    with tabs[3]:
        render_network()

