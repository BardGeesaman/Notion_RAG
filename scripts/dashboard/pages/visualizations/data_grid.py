from __future__ import annotations

import pandas as pd
import streamlit as st

from scripts.dashboard.db_session import db_session
from amprenta_rag.database.models import Dataset, Feature

try:
    # Provided by the `streamlit-aggrid` package
    from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode  # type: ignore
except Exception:  # pragma: no cover
    AgGrid = None  # type: ignore[assignment]
    GridOptionsBuilder = None  # type: ignore[assignment]
    GridUpdateMode = None  # type: ignore[assignment]


def _load_datasets(limit: int = 200):
    with db_session() as db:
        rows = (
            db.query(Dataset.id, Dataset.name, Dataset.omics_type, Dataset.description, Dataset.created_at)
            .order_by(Dataset.created_at.desc())
            .limit(limit)
            .all()
        )
        data = []
        for r in rows:
            data.append(
                {
                    "id": str(r.id),
                    "name": r.name,
                    "omics_type": r.omics_type,
                    "description": r.description or "",
                    "created_at": r.created_at.isoformat() if r.created_at else "",
                }
            )
        return pd.DataFrame(data)


def _load_features(limit: int = 500):
    with db_session() as db:
        rows = (
            db.query(Feature.id, Feature.name, Feature.feature_type, Feature.normalized_name, Feature.created_at)
            .order_by(Feature.created_at.desc())
            .limit(limit)
            .all()
        )
        data = []
        for r in rows:
            data.append(
                {
                    "id": str(r.id),
                    "name": r.name,
                    "feature_type": r.feature_type,
                    "normalized_name": r.normalized_name or "",
                    "created_at": r.created_at.isoformat() if r.created_at else "",
                }
            )
        return pd.DataFrame(data)


def render() -> None:
    st.header("AG-Grid Data Browser")
    st.caption("Browse and edit entities using interactive AG-Grid tables.")

    if AgGrid is None or GridOptionsBuilder is None or GridUpdateMode is None:
        st.error("❌ Data Grid requires `streamlit-aggrid` (module `st_aggrid`) which is not installed.")
        st.code("pip install streamlit-aggrid", language="bash")
        st.info("Showing a basic (non-editable) fallback table instead.")

    entity = st.selectbox("Entity", ["Datasets", "Features"], index=0)
    page_size = st.slider("Page size", min_value=20, max_value=200, value=50, step=10)

    if entity == "Datasets":
        df = _load_datasets(limit=2000)
        editable_cols = ["description"]
    else:
        df = _load_features(limit=5000)
        editable_cols = ["normalized_name"]

    if df.empty:
        st.info("No records found.")
        return

    if AgGrid is None or GridOptionsBuilder is None or GridUpdateMode is None:
        st.dataframe(df, use_container_width=True, hide_index=True)
        csv = df.to_csv(index=False)
        st.download_button(
            label="📥 Export CSV",
            data=csv,
            file_name=f"{entity.lower().replace(' ', '_')}.csv",
            mime="text/csv",
        )
        return

    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=page_size)
    gb.configure_default_column(editable=False, filter=True, sortable=True, resizable=True)
    for col in editable_cols:
        gb.configure_column(col, editable=True)
    grid_options = gb.build()

    st.subheader(f"{entity} ({len(df)} rows)")
    grid_response = AgGrid(
        df,
        gridOptions=grid_options,
        update_mode=GridUpdateMode.VALUE_CHANGED,
        enable_enterprise_modules=False,
        fit_columns_on_grid_load=True,
        height=600,
    )

    if grid_response and grid_response.get("data") is not None:
        edited_df = pd.DataFrame(grid_response["data"])
        if not edited_df.equals(df):
            st.warning("Edits are local only (not persisted).")

    csv = df.to_csv(index=False)
    st.download_button(
        label="📥 Export CSV",
        data=csv,
        file_name=f"{entity.lower().replace(' ', '_')}.csv",
        mime="text/csv",
    )

