"""Experiment Optimizer page (Bayesian optimization recommendations)."""

from __future__ import annotations

import json
from typing import Any, Dict, List

import streamlit as st

from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds


DEFAULT_TESTED = [
    {"features": [0.0, 0.0], "activity": 0.1},
    {"features": [0.0, 1.0], "activity": 0.2},
    {"features": [1.0, 0.0], "activity": 0.25},
    {"features": [1.0, 1.0], "activity": 0.6},
    {"features": [2.0, 2.0], "activity": 0.9},
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


def render_experiment_optimizer_page() -> None:
    """Render the Experiment Optimizer page."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧠 Experiment Optimizer")
    st.caption("BoTorch-based Bayesian optimization to recommend next compounds to test.")

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
            recs = recommend_next_compounds(
                tested_compounds=tested,
                candidate_pool=pool,
                batch_size=int(batch_size),
                objectives=["potency"],
            )
            st.session_state["optimizer_recs"] = recs
        except ImportError as e:
            st.error(str(e))
            st.info("Install dependencies and restart: botorch, ax-platform, torch.")
            return
        except Exception as e:  # noqa: BLE001
            st.error(f"Recommendation failed: {e}")
            return

    recs = st.session_state.get("optimizer_recs") or []
    st.markdown("### Recommendations")
    if not recs:
        st.info("Run Recommend to see ranked suggestions.")
        return

    st.dataframe(recs, use_container_width=True, hide_index=True)


__all__ = ["render_experiment_optimizer_page"]


