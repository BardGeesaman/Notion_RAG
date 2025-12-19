import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset, Program
from scripts.dashboard.auth import require_auth
from scripts.dashboard.db_session import db_session


def render_coverage_page():
    require_auth()
    st.header("Multi-Omics Coverage Map")
    st.info(
        """Dark cells = more datasets; empty = data gap. Status now shown by explicit ingestion_status, not inferred."""
    )
    with db_session() as db:
        datasets = db.query(Dataset).all()
        {p.id: p.name for p in db.query(Program).all()}
        df = pd.DataFrame(
            [
                {
                    "program": next((p.name for p in getattr(d, "programs", [])[:1]), None),
                    "disease": d.disease[0] if getattr(d, "disease", None) else None,
                    "omics": d.omics_type,
                    "status": getattr(d, "ingestion_status", "pending"),
                }
                for d in datasets
            ]
        )
        # Handle empty DataFrame gracefully
        if df.empty:
            st.warning("No datasets found. Ingest some data to see the coverage map.")
            return

        filter_status = st.selectbox("Filter by ingestion status", ["All"] + sorted(df["status"].dropna().unique().tolist()))
        filter_prog = st.selectbox("Filter by program", ["All"] + sorted(df["program"].dropna().unique().tolist()))
        filter_dis = st.selectbox("Filter by disease", ["All"] + sorted(df["disease"].dropna().unique().tolist()))
        abs_toggle = st.toggle("Show absolute counts", value=True)
        view = df
        if filter_status != "All":
            view = view[view["status"] == filter_status]
        if filter_prog != "All":
            view = view[view["program"] == filter_prog]
        if filter_dis != "All":
            view = view[view["disease"] == filter_dis]
        pivot = view.pivot_table(index="program", columns="omics", aggfunc="size", fill_value=0)
        if not abs_toggle:
            pivot = pivot.map(lambda x: 1 if x > 0 else 0)
        st.subheader("Program × Omics")
        st.dataframe(pivot, width='stretch')
        st.subheader("Disease × Omics")
        pivot2 = view.pivot_table(index="disease", columns="omics", aggfunc="size", fill_value=0)
        if not abs_toggle:
            pivot2 = pivot2.map(lambda x: 1 if x > 0 else 0)
        st.dataframe(pivot2, width='stretch')
        st.caption("Status is now explicit (pending, in_progress, complete, failed)")
