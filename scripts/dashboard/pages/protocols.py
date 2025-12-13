"""Protocol Management - Create and manage experimental protocols."""
import streamlit as st
import json
from datetime import datetime
from amprenta_rag.database.models import Protocol, ExperimentProtocol, Experiment
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.session import get_current_user


def render_protocols_page():
    st.title("📋 Protocol Management")

    tab1, tab2, tab3 = st.tabs(["Browse Protocols", "Create Protocol", "Link to Experiment"])

    with tab1:
        render_browse_tab()

    with tab2:
        render_create_tab()

    with tab3:
        render_link_tab()


def render_browse_tab():
    st.subheader("Protocol Library")

    with db_session() as db:
        protocols = db.query(Protocol).filter(Protocol.is_active == True).order_by(Protocol.name).all()

        if not protocols:
            st.info("No protocols found. Create your first protocol!")
            return

        for p in protocols:
            with st.expander(f"**{p.name}** (v{p.version}) - {p.category or 'Uncategorized'}"):
                st.markdown(p.description or "_No description_")

                if p.steps:
                    st.markdown("**Steps:**")
                    for step in p.steps:
                        st.markdown(f"  {step.get('order', '?')}. {step.get('title', 'Untitled')}")

                if p.materials:
                    st.markdown("**Materials:**")
                    for mat in p.materials:
                        st.markdown(f"  - {mat.get('name')}: {mat.get('quantity', '')} {mat.get('unit', '')}")

                col1, col2 = st.columns(2)
                with col1:
                    if st.button("📝 New Version", key=f"version_{p.id}"):
                        st.session_state[f"clone_{p.id}"] = True
                        st.rerun()


def render_create_tab():
    st.subheader("Create New Protocol")
    user = get_current_user()

    with st.form("create_protocol"):
        name = st.text_input("Protocol Name*")
        category = st.selectbox("Category", ["extraction", "assay", "analysis", "culture", "purification", "other"])
        description = st.text_area("Description")

        st.markdown("**Steps** (JSON array)")
        steps_json = st.text_area("Steps", value='[{"order": 1, "title": "Step 1", "instructions": "", "duration": ""}]', height=100)

        st.markdown("**Materials** (JSON array)")
        materials_json = st.text_area("Materials", value='[{"name": "", "quantity": "", "unit": ""}]', height=100)

        submitted = st.form_submit_button("Create Protocol", type="primary")

        if submitted:
            if not name:
                st.error("Protocol name is required")
                return

            try:
                steps = json.loads(steps_json) if steps_json.strip() else None
                materials = json.loads(materials_json) if materials_json.strip() else None
            except json.JSONDecodeError as e:
                st.error(f"Invalid JSON: {e}")
                return

            with db_session() as db:
                protocol = Protocol(
                    name=name,
                    category=category,
                    description=description,
                    steps=steps,
                    materials=materials,
                    created_by_id=user.get("id") if user and user.get("id") != "test" else None
                )
                db.add(protocol)
                db.commit()
                st.success(f"Protocol '{name}' created successfully!")


def render_link_tab():
    st.subheader("Link Protocol to Experiment")

    with db_session() as db:
        protocols = db.query(Protocol).filter(Protocol.is_active == True).all()
        experiments = db.query(Experiment).order_by(Experiment.name).limit(100).all()

        if not protocols:
            st.warning("No protocols available")
            return
        if not experiments:
            st.warning("No experiments available")
            return

        protocol_options = {f"{p.name} (v{p.version})": p.id for p in protocols}
        experiment_options = {e.name: e.id for e in experiments}

        with st.form("link_protocol"):
            selected_exp = st.selectbox("Experiment", list(experiment_options.keys()))
            selected_proto = st.selectbox("Protocol", list(protocol_options.keys()))
            notes = st.text_area("Execution Notes")
            executed = st.checkbox("Mark as executed")

            if st.form_submit_button("Link Protocol"):
                link = ExperimentProtocol(
                    experiment_id=experiment_options[selected_exp],
                    protocol_id=protocol_options[selected_proto],
                    notes=notes if notes else None,
                    executed_at=datetime.utcnow() if executed else None
                )
                db.add(link)
                db.commit()
                st.success("Protocol linked to experiment!")


if __name__ == "__main__":
    render_protocols_page()
