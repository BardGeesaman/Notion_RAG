import streamlit as st

from amprenta_rag.database.models import Dataset


def render_dataset_ingest_summary(dataset: Dataset):
    if not dataset:
        st.success("✅ Dataset ingested successfully (but could not fetch details).")
        return
    st.success("✅ Ingestion Complete!")
    st.markdown(f"**Dataset:** {dataset.name if hasattr(dataset, 'name') else ''}")
    st.markdown(f"**ID:** `{dataset.id}`")
    st.markdown(f"**Omics Type:** {getattr(dataset, 'omics_type', '')}")
    prog = getattr(dataset, "programs", [None])[0]
    exp = getattr(dataset, "experiments", [None])[0]
    st.markdown(f"**Program:** {prog.name if prog else '(None)'}")
    st.markdown(f"**Experiment:** {exp.name if exp else '(None)'}")
    st.markdown(f"**Ingestion Status:** {getattr(dataset, 'ingestion_status', 'N/A')}")

    st.markdown("---")
    st.markdown("### 🎯 Next Steps")
    st.info(
        f"""
    **View your dataset:**
    - Go to **Datasets** page in the sidebar
    - Search for dataset ID: `{dataset.id}`

    **Explore the data:**
    - **Features** page - View extracted features
    - **Signatures** page - Match against known signatures
    - **Cross-Omics** page - Cross-reference with other datasets
    """
    )
