import streamlit as st

from amprenta_rag.database.models import Dataset
from amprenta_rag.reporting.evidence_report import generate_dataset_report
from scripts.dashboard.db_session import db_session


def dataset_selector():
    with db_session() as db:
        dsets = db.query(Dataset).order_by(Dataset.created_at.desc()).all()
        # Extract data inside session before it closes
        id_names = [(str(d.id), f"{d.name} ({d.omics_type})") for d in dsets]

    ds_id = st.selectbox(
        "Select a dataset", [id_ for id_, _ in id_names], format_func=lambda x: dict(id_names).get(x, x)
    )
    return ds_id


def render_evaluation_wizard():
    st.title("🧭 Dataset Evaluation Wizard")
    if "wizard_dataset_id" not in st.session_state:
        st.session_state["wizard_dataset_id"] = None
    if "current_step" not in st.session_state:
        st.session_state["current_step"] = 1
    step = st.session_state["current_step"]
    N_STEPS = 5
    st.progress((step - 1) / N_STEPS, text=f"Step {step} of {N_STEPS}")
    # Step navigation
    nav_cols = st.columns(N_STEPS)
    for idx in range(N_STEPS):
        if nav_cols[idx].button(f"{idx+1}"):
            st.session_state["current_step"] = idx + 1
            st.rerun()
    # Step logic
    if step == 1:
        st.header("Step 1: Select a Dataset")
        ds_id = st.session_state["wizard_dataset_id"] or dataset_selector()
        st.session_state["wizard_dataset_id"] = ds_id
        with db_session() as db:
            ds = db.query(Dataset).filter(Dataset.id == ds_id).first()
            if ds:
                # Extract data inside session
                ds_name = ds.name
                ds_omics = ds.omics_type
                ds_program = ds.programs[0].name if ds.programs else "(None)"
                ds_status = ds.ingestion_status

        if ds_id:
            st.markdown(f"**Name:** {ds_name}")
            st.markdown(f"**Omics:** {ds_omics}")
            st.markdown(f"**Program:** {ds_program}")
            st.markdown(f"**Ingestion Status:** {ds_status}")
        if st.button("Next: QC & Coverage →", type="primary"):
            st.session_state["current_step"] = 2
            st.rerun()
    elif step == 2:
        st.header("Step 2: Data Review & Coverage")
        ds_id = st.session_state["wizard_dataset_id"]
        with db_session() as db:
            ds = db.query(Dataset).filter(Dataset.id == ds_id).first()
            if ds:
                # Extract data inside session
                ds_ingestion_status = ds.ingestion_status
                ds_name = ds.name
                feature_count = len(ds.features)
            else:
                ds_ingestion_status = None
                ds_name = None
                feature_count = 0

        if ds_ingestion_status:
            st.markdown(f"##### Ingestion Status: {ds_ingestion_status}")
            st.markdown(f"**Dataset:** {ds_name}")
            st.markdown(f"**Features:** {feature_count}")
            st.info("💡 Use the **Datasets** or **Data Management** page in the sidebar to view full details")
        if st.button("Next: Discovery/Signatures →", type="primary"):
            st.session_state["current_step"] = 3
            st.rerun()
        if st.button("← Back", type="secondary"):
            st.session_state["current_step"] = 1
            st.rerun()
    elif step == 3:
        st.header("Step 3: Run Discovery & Signature Matching")
        ds_id = st.session_state["wizard_dataset_id"]
        st.info(f"💡 Go to **Signatures** page in the sidebar to discover signatures for dataset: `{ds_id}`")
        st.markdown("[This will let you assess and link signatures against your dataset.]")
        if st.button("Next: Generate Evidence Report →", type="primary"):
            st.session_state["current_step"] = 4
            st.rerun()
        if st.button("← Back", type="secondary"):
            st.session_state["current_step"] = 2
            st.rerun()
    elif step == 4:
        st.header("Step 4: Generate Evidence Report")
        ds_id = st.session_state["wizard_dataset_id"]
        if st.button("Create Evidence Report Summary", key="gen_report"):
            with st.spinner("Building report..."):
                try:
                    report = generate_dataset_report(ds_id)
                    st.session_state["last_report"] = report
                except Exception as e:
                    st.error(str(e))
        if "last_report" in st.session_state:
            report = st.session_state["last_report"]
            st.markdown(f"**Summary for {report.entity_id}:**")
            for s in report.sections:
                st.markdown(f"**{s.title}**\n\n{s.summary_text}")
        st.info(f"💡 Full evidence report for dataset: `{ds_id}`")
        if st.button("Next: Ask the Assistant →", type="primary"):
            st.session_state["current_step"] = 5
            st.rerun()
        if st.button("← Back", type="secondary"):
            st.session_state["current_step"] = 3
            st.rerun()
    elif step == 5:
        st.header("Step 5: Ask the Assistant about this Dataset")
        ds_id = st.session_state["wizard_dataset_id"]
        st.info("💡 Go to the **Chat** page in the sidebar to ask questions about this dataset")
        st.markdown(
            f"""
        **Example questions:**
        - What are the key findings for dataset `{ds_id}`?
        - How does this dataset relate to ALS?
        - What signatures are relevant to this dataset?
        """
        )
        if st.button("Restart Wizard", type="secondary"):
            st.session_state["current_step"] = 1
            st.session_state["wizard_dataset_id"] = None
            st.rerun()
