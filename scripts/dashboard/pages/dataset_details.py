"""
Dataset Details page.
"""

from __future__ import annotations

import json
import os
from typing import Any, Optional
from uuid import UUID

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


def patch_dataset(dataset_id: str, payload: dict) -> httpx.Response:
    """Patch a dataset via the API."""
    with httpx.Client(timeout=15) as client:
        return client.patch(f"{DATASET_ENDPOINT}/{dataset_id}", json=payload)


def _parse_csv_list(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


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

    # Paginated feature display
    page_size = st.session_state.get("feature_page_size", 10)
    sample_rows = []
    for ftype, ids in feature_ids.items():
        for fid in ids[:page_size]:
            sample_rows.append({"Feature Type": ftype, "Feature ID": fid})
    if sample_rows:
        st.table(pd.DataFrame(sample_rows))

        total_features = sum(len(ids) for ids in feature_ids.values())
        if total_features > page_size:
            if st.button(f"Show more (viewing {page_size} of {total_features})"):
                st.session_state["feature_page_size"] = page_size + 20
                st.rerun()


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
    version = int(dataset.get("version") or 1)
    st.session_state["dataset_version"] = version

    st.header(dataset.get("name") or "Dataset")
    st.caption(f"v{version}")
    col1, col2 = st.columns([2, 1])
    with col1:
        render_metadata(dataset)
    with col2:
        render_feature_summary(dataset)

    st.markdown("---")
    st.subheader("Edit Dataset")

    with st.form("dataset_update_form"):
        name = st.text_input("Name", value=dataset.get("name") or "")
        description = st.text_area("Description", value=dataset.get("description") or "")
        organism_csv = st.text_input(
            "Organism (comma-separated)",
            value=", ".join(dataset.get("organism") or []),
        )
        sample_type_csv = st.text_input(
            "Sample type (comma-separated)",
            value=", ".join(dataset.get("sample_type") or []),
        )
        disease_csv = st.text_input(
            "Disease (comma-separated)",
            value=", ".join(dataset.get("disease") or []),
        )
        submitted = st.form_submit_button("Save", type="primary")

    if submitted:
        payload: dict[str, Any] = {
            "version": int(st.session_state.get("dataset_version", version)),
            "name": name.strip() or None,
            "description": description.strip() or None,
            "organism": _parse_csv_list(organism_csv) if organism_csv.strip() else [],
            "sample_type": _parse_csv_list(sample_type_csv) if sample_type_csv.strip() else [],
            "disease": _parse_csv_list(disease_csv) if disease_csv.strip() else [],
        }
        payload = {k: v for k, v in payload.items() if v is not None}

        try:
            resp = patch_dataset(str(dataset.get("id") or ""), payload)
        except httpx.HTTPError as e:  # noqa: BLE001
            st.error(f"Failed to update dataset: {e}")
            return

        if resp.status_code == 409:
            st.error("This dataset was modified. Please reload.")
            if st.button("🔄 Reload", type="secondary"):
                fetch_dataset.clear()
                st.rerun()
            return

        if resp.status_code == 404:
            st.error("Dataset not found.")
            return

        if resp.status_code >= 500:
            st.error("Server error - please try again.")
            return

        try:
            resp.raise_for_status()
        except httpx.HTTPError as e:  # noqa: BLE001
            st.error(f"Update failed: {e}")
            return

        updated = resp.json()
        st.session_state["dataset_version"] = int(updated.get("version") or (version + 1))
        st.success("Saved.")
        fetch_dataset.clear()
        st.rerun()


def render_selection():
    initial_id = st.query_params.get("dataset_id")
    dataset_id = st.text_input("Dataset ID", value=initial_id or "")
    load_clicked = st.button("Load Dataset", type="primary")
    if dataset_id:
        try:
            UUID(dataset_id)
        except ValueError:
            st.error("Invalid UUID format. Please enter a valid dataset ID.")
            return None, False
    return dataset_id, load_clicked or bool(initial_id)


def render_page():
    st.set_page_config(page_title="Dataset Details", layout="wide")
    # Breadcrumb navigation
    cols = st.columns([1, 1, 8])
    with cols[0]:
        if st.button("Catalog"):
            st.switch_page("pages/external_catalog.py")
    with cols[1]:
        st.markdown(" → **Dataset**")

    st.title("Dataset Details")

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

