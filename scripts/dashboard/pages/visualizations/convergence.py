from __future__ import annotations

from collections import defaultdict
from typing import Dict, List

import pandas as pd
import plotly.express as px
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Feature, Signature
from scripts.dashboard.db_session import db_session


def _gather_counts(db: Session) -> pd.DataFrame:
    signatures: List[Signature] = db.query(Signature).all()
    if not signatures:
        return pd.DataFrame()

    counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for sig in signatures:
        sig_label = sig.name or str(sig.id)
        for feat in sig.features:
            omics = getattr(feat, "omics_type", None) or _infer_omics(feat)
            counts[sig_label][omics] += 1

    records = []
    for sig_label, omics_counts in counts.items():
        for omics_type, count in omics_counts.items():
            records.append({"signature": sig_label, "omics": omics_type, "count": count})

    return pd.DataFrame(records)


def _infer_omics(feat: Feature) -> str:
    name = (feat.name or "").lower()
    # Rough heuristics
    if any(k in name for k in ["hsa", "gene", "ensg", "tp53", "akt"]):
        return "transcriptomics"
    if any(k in name for k in ["P0", "Q9", "uniprot", "prot"]):
        return "proteomics"
    if any(k in name for k in ["met", "gly", "acid", "glc", "lactate"]):
        return "metabolomics"
    if any(k in name for k in ["cer", "lipid", "sm", "pc(", "pe("]):
        return "lipidomics"
    return "unknown"


def render() -> None:
    st.header("Cross-Omics Convergence")
    st.caption("Feature counts per omics type across signatures.")

    with db_session() as db:
        df = _gather_counts(db)

    if df.empty:
        st.info("No signature features available.")
        return

    fig = px.bar(
        df,
        x="signature",
        y="count",
        color="omics",
        barmode="group",
        title="Omics composition per signature",
    )
    fig.update_layout(xaxis_tickangle=45)
    st.plotly_chart(fig, use_container_width=True)

    st.download_button(
        "Download Counts (CSV)",
        data=df.to_csv(index=False),
        file_name="cross_omics_convergence.csv",
        mime="text/csv",
    )


__all__ = ["render"]
