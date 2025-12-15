"""MOA Inference dashboard page."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.analysis.moa_inference import infer_moa
from amprenta_rag.database.models import Compound, Dataset
from amprenta_rag.database.session import db_session


def render_moa_inference_page() -> None:
    st.header("🧭 MOA Inference")
    st.caption("Infer mechanism of action by combining HTS bioactivity, multi-omics concordance, and pathway signals.")

    with db_session() as db:
        compounds = db.query(Compound).order_by(Compound.created_at.desc()).limit(200).all()
        datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).limit(200).all()

    if not compounds:
        st.info("No compounds available.")
        return

    comp_map = {c.compound_id: c for c in compounds}
    selected_comp = st.selectbox("Select compound", list(comp_map.keys()))
    ds_map = {d.name: d for d in datasets}
    selected_ds = st.multiselect("Select datasets (optional)", list(ds_map.keys()))
    dataset_ids = [ds_map[n].id for n in selected_ds]

    if st.button("Run MOA Inference", type="primary"):
        try:
            candidates = infer_moa(comp_map[selected_comp].id, dataset_ids)
        except Exception as exc:  # pragma: no cover - UI fallback
            st.error(f"Failed to infer MOA: {exc}")
            return

        if not candidates:
            st.info("No MOA candidates generated.")
            return

        rows = []
        for c in candidates:
            rows.append(
                {
                    "Rank": c.rank,
                    "Candidate": c.candidate_id,
                    "Type": c.type,
                    "Probability": c.probability,
                }
            )
        df = pd.DataFrame(rows)
        st.subheader("Ranked Candidates")
        st.dataframe(df, use_container_width=True, hide_index=True)

        st.markdown("---")
        st.subheader("Evidence Breakdown")
        for c in candidates:
            with st.expander(f"{c.rank}. {c.candidate_id} (p={c.probability:.2f})"):
                evid_rows = [
                    {
                        "Feature": e.feature_name,
                        "Value": e.value,
                        "Weight": e.weight,
                    }
                    for e in c.contributions
                ]
                st.table(pd.DataFrame(evid_rows))

