"""
Dataset Details page.
"""

from __future__ import annotations

import json
import os
from typing import Optional

import httpx
import pandas as pd
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
DATASET_ENDPOINT = f"{API_BASE}/api/v1/datasets"


@st.cache_data(ttl=30)
def fetch_dataset(dataset_id: str) -> Optional[dict]:
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(f"{DATASET_ENDPOINT}/{dataset_id}")
            if resp.status_code == 404:
                st.error("Dataset not found.")
                return None
            if resp.status_code >= 500:
                st.error("Server error - please try again.")
                return None
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load dataset: {e}")
        return None


def get_external_link(source: Optional[str], accession: Optional[str]) -> Optional[str]:
    if not source or not accession:
        return None
    src = source.lower()
    if src == "geo":
        return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    if src == "pride":
        return f"https://www.ebi.ac.uk/pride/archive/projects/{accession}"
    if src == "metabolights":
        return f"https://www.ebi.ac.uk/metabolights/{accession}"
    if src == "mw":
        return f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={accession}"
    return None


def render_feature_summary(dataset: dict):
    feature_ids = dataset.get("feature_ids") or {}
    if not feature_ids:
        st.info("No feature summary available.")
        return

    # Count by type
    counts = {ftype: len(ids) for ftype, ids in feature_ids.items()}
    count_df = pd.DataFrame({"Feature Type": list(counts.keys()), "Count": list(counts.values())})
    st.subheader("Feature Summary")
    st.bar_chart(count_df.set_index("Feature Type"))

    # Sample values
    sample_rows = []
    for ftype, ids in feature_ids.items():
        for fid in ids[:10]:
            sample_rows.append({"Feature Type": ftype, "Feature ID": fid})
    if sample_rows:
        st.table(pd.DataFrame(sample_rows))


def render_metadata(dataset: dict):
    ext = dataset.get("external_ids") or {}
    accession = ext.get("accession") or ext.get("study_id") or ext.get("id")
    source = dataset.get("data_origin") or dataset.get("source")

    st.subheader("Metadata")
    st.write(f"**Source:** {source or 'Unknown'}")
    st.write(f"**Accession:** {accession or 'N/A'}")
    st.write(f"**Import Date:** {dataset.get('created_at')}")
    st.write(f"**Status:** {dataset.get('ingestion_status') or 'unknown'}")
    if dataset.get("description"):
        st.markdown(f"**Description:** {dataset['description']}")

    ext_link = get_external_link(source, accession)
    if ext_link:
        st.link_button("Open in Repository", ext_link, use_container_width=True)

    # Export metadata
    meta_json = json.dumps(dataset, indent=2, default=str)
    st.download_button("Export metadata (JSON)", meta_json, file_name="dataset_metadata.json", mime="application/json")


def render_dataset_details(dataset: dict):
    st.header(dataset.get("name") or "Dataset")
    col1, col2 = st.columns([2, 1])
    with col1:
        render_metadata(dataset)
    with col2:
        render_feature_summary(dataset)


def render_selection():
    initial_id = st.query_params.get("dataset_id")
    dataset_id = st.text_input("Dataset ID", value=initial_id or "")
    load_clicked = st.button("Load Dataset", type="primary")
    return dataset_id, load_clicked or bool(initial_id)


def render_page():
    st.set_page_config(page_title="Dataset Details", layout="wide")
    st.title("Dataset Details")
    # Back link to catalog
    if st.button("← Back to Catalog"):
        st.switch_page("pages/external_catalog.py")

    dataset_id, should_load = render_selection()
    if not should_load or not dataset_id:
        st.info("Enter a dataset ID to view details.")
        return

    with st.spinner("Loading dataset..."):
        dataset = fetch_dataset(dataset_id)
    if not dataset:
        return

    render_dataset_details(dataset)


if __name__ == "__main__":
    render_page()

