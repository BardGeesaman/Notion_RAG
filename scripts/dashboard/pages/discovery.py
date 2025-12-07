import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from amprenta_rag.domain.signatures_discovery import DiscoveryDatasetSummary
from amprenta_rag.signatures.signature_discovery import discover_signatures_from_datasets
from scripts.dashboard.db_session import db_session


def render_discovery_page():
    st.header("Signature Discovery")
    with db_session() as db:
        # Filters
        omics_types = sorted(set(d.omics_type for d in db.query(Dataset).all() if d.omics_type))
        omics_type = st.selectbox("Omics type", omics_types)
        disease = st.text_input("Disease filter (optional)", value="")
        matrix = st.text_input("Matrix filter (optional)", value="")
        min_support = st.number_input("Minimum support (datasets)", min_value=1, value=2)
        min_overlap = st.slider("Minimum overlap", min_value=0.0, max_value=1.0, value=0.3, step=0.05)
        if st.button("Run discovery"):
            # Construct summaries:
            ds_query = db.query(Dataset).filter(Dataset.omics_type == omics_type)
            if disease:
                ds_query = ds_query.filter(Dataset.disease.contains([disease]))
            if matrix:
                ds_query = ds_query.filter(Dataset.matrix == matrix)
            datasets = ds_query.all()
            dsums = [
                DiscoveryDatasetSummary(
                    dataset_id=d.id,
                    omics_type=d.omics_type,
                    disease=d.disease[0] if d.disease else None,
                    matrix=d.matrix,
                    features={f.name for f in d.features} if d.features else set(),
                    directions=None,
                )
                for d in datasets
            ]
            sigs = discover_signatures_from_datasets(dsums, min_support=min_support, min_overlap=min_overlap)
            table = [
                {
                    "Name": sig.name,
                    "Components": len(sig.components),
                    "Support": sig.support,
                    "Provenance": ", ".join(sig.provenance.get("dataset_ids", [])),
                    "Description": "Auto-discovered signature",
                }
                for sig in sigs
            ]
            df = pd.DataFrame(table)
            st.dataframe(df, use_container_width=True, hide_index=True)
            for i, sig in enumerate(sigs):
                # Download TSV
                if st.download_button(
                    f"Download {sig.name} TSV",
                    data="\n".join([c.feature for c in sig.components]),
                    file_name=f"{sig.name}.tsv",
                    key=f"tsv_{i}",
                ):
                    pass
                # Evidence report link (simulate or real route)
                if st.button(f"Open in Evidence Report: {sig.name}", key=f"rep_{i}"):
                    # Set session variable or navigate, depending on your dashboard router.
                    st.session_state["selected_signature_id"] = str(sig.name)
