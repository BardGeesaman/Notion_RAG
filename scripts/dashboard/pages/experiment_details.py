"""
Experiment Details page (API-backed).

This page loads an experiment from the FastAPI service and allows editing core
fields. Updates use optimistic locking by sending the current `version` and
handling HTTP 409 conflicts with a reload prompt.
"""

from __future__ import annotations

import os
from typing import Any, Optional
from uuid import UUID

import httpx
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
EXPERIMENT_ENDPOINT = f"{API_BASE}/api/v1/experiments"


@st.cache_data(ttl=30)
def fetch_experiment(experiment_id: str) -> Optional[dict]:
    """Fetch an experiment from the API."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(f"{EXPERIMENT_ENDPOINT}/{experiment_id}")
            if resp.status_code == 404:
                st.error("Experiment not found.")
                return None
            if resp.status_code >= 500:
                st.error("Server error - please try again.")
                return None
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load experiment: {e}")
        return None


def patch_experiment(experiment_id: str, payload: dict) -> httpx.Response:
    """Patch an experiment via the API."""
    with httpx.Client(timeout=15) as client:
        return client.patch(f"{EXPERIMENT_ENDPOINT}/{experiment_id}", json=payload)


def _parse_csv_list(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


def render_page() -> None:
    st.set_page_config(page_title="Experiment Details", layout="wide")
    st.title("Experiment Details")

    initial_id = st.query_params.get("experiment_id")
    experiment_id = st.text_input("Experiment ID", value=initial_id or "")
    load_clicked = st.button("Load Experiment", type="primary")

    if experiment_id:
        try:
            UUID(experiment_id)
        except ValueError:
            st.error("Invalid UUID format. Please enter a valid experiment ID.")
            return

    should_load = load_clicked or bool(initial_id)
    if not should_load or not experiment_id:
        st.info("Enter an experiment ID to view details.")
        return

    with st.spinner("Loading experiment..."):
        exp = fetch_experiment(experiment_id)
    if not exp:
        return

    version = int(exp.get("version") or 1)
    st.session_state["experiment_version"] = version

    header = exp.get("name") or "Experiment"
    st.header(header)
    st.caption(f"v{version}")

    with st.form("experiment_update_form"):
        name = st.text_input("Name", value=exp.get("name") or "")
        exp_type = st.text_input("Type", value=exp.get("type") or "")
        description = st.text_area("Description", value=exp.get("description") or "")

        disease_csv = st.text_input(
            "Disease (comma-separated)",
            value=", ".join(exp.get("disease") or []),
        )
        matrix_csv = st.text_input(
            "Matrix (comma-separated)",
            value=", ".join(exp.get("matrix") or []),
        )
        model_systems_csv = st.text_input(
            "Model systems (comma-separated)",
            value=", ".join(exp.get("model_systems") or []),
        )

        submitted = st.form_submit_button("Save", type="primary")

    if submitted:
        payload: dict[str, Any] = {
            "version": int(st.session_state.get("experiment_version", version)),
            "name": name.strip() or None,
            "type": exp_type.strip() or None,
            "description": description.strip() or None,
            "disease": _parse_csv_list(disease_csv) if disease_csv.strip() else [],
            "matrix": _parse_csv_list(matrix_csv) if matrix_csv.strip() else [],
            "model_systems": _parse_csv_list(model_systems_csv) if model_systems_csv.strip() else [],
        }
        # Remove explicit None fields to avoid overwriting unless user intended.
        payload = {k: v for k, v in payload.items() if v is not None}

        try:
            resp = patch_experiment(experiment_id, payload)
        except httpx.HTTPError as e:  # noqa: BLE001
            st.error(f"Failed to update experiment: {e}")
            return

        if resp.status_code == 409:
            st.error("This experiment was modified. Please reload.")
            if st.button("🔄 Reload", type="secondary"):
                fetch_experiment.clear()
                st.rerun()
            return

        if resp.status_code == 404:
            st.error("Experiment not found.")
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
        new_version = int(updated.get("version") or (version + 1))
        st.session_state["experiment_version"] = new_version
        st.success("Saved.")
        fetch_experiment.clear()
        st.rerun()


if __name__ == "__main__":
    render_page()


