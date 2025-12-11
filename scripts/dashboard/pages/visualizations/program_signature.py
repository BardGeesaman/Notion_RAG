from __future__ import annotations

from typing import Dict, Set

import pandas as pd
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Program
from scripts.dashboard.db_session import db_session
from scripts.dashboard.utils.viz_helpers import create_heatmap


def _build_matrix(db: Session) -> pd.DataFrame:
    programs = db.query(Program).all()
    if not programs:
        return pd.DataFrame()

    signature_names: Set[str] = set()
    data: Dict[str, Dict[str, int]] = {}

    for prog in programs:
        prog_label = prog.name or str(prog.id)
        row: Dict[str, int] = {}
        for sig in prog.signatures:
            sig_label = sig.name or str(sig.id)
            signature_names.add(sig_label)
            row[sig_label] = 1
        data[prog_label] = row

    if not signature_names:
        return pd.DataFrame()

    df = pd.DataFrame(0, index=sorted(data.keys()), columns=sorted(signature_names))
    for prog_label, row in data.items():
        for sig_label in row:
            df.at[prog_label, sig_label] = 1
    return df


def render() -> None:
    st.header("Program-Signature Map")
    st.caption("Matrix of program to signature links (1 = linked, 0 = not linked).")

    with db_session() as db:
        mat = _build_matrix(db)

    if mat.empty:
        st.info("No program-signature relationships found.")
        return

    fig = create_heatmap(
        mat,
        title="Program × Signature",
        colorscale="Blues",
        cluster=False,
    )
    st.plotly_chart(fig, width='stretch')

    st.download_button(
        "Download Matrix (CSV)",
        data=mat.to_csv(),
        file_name="program_signature_matrix.csv",
        mime="text/csv",
    )


__all__ = ["render"]
