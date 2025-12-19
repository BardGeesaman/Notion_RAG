"""Candidate Selection page for TPP-based compound scoring."""
from __future__ import annotations

from uuid import UUID

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Compound, TargetProductProfile, CandidateNomination
from amprenta_rag.chemistry.candidate_selection import score_compound_against_tpp, nominate_candidate
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_candidate_selection_page() -> None:
    """Render the Candidate Selection page."""
    st.header("🎯 Candidate Selection")
    st.markdown("Define Target Product Profiles (TPP) and score compounds against them.")

    tab1, tab2, tab3 = st.tabs(["Define TPP", "Score Compound", "Nominations"])

    with tab1:
        render_tpp_tab()
    with tab2:
        render_score_tab()
    with tab3:
        render_nominations_tab()


def render_tpp_tab() -> None:
    """Render the Define TPP tab."""
    st.subheader("Target Product Profiles")

    user = get_current_user()

    with db_session() as db:
        # List existing TPPs
        tpps = db.query(TargetProductProfile).order_by(TargetProductProfile.created_at.desc()).all()

        if tpps:
            st.markdown("### Existing TPPs")
            for tpp in tpps:
                with st.expander(f"**{tpp.name}** {'✅' if tpp.is_active else '❌'}"):
                    st.markdown(f"**Description:** {tpp.description or 'No description'}")
                    st.markdown(f"**Criteria:** {len(tpp.criteria) if tpp.criteria else 0} criteria")
                    if tpp.criteria:
                        st.json(tpp.criteria)

        st.markdown("---")
        st.markdown("### Create New TPP")

        with st.form("create_tpp"):
            name = st.text_input("TPP Name*", placeholder="e.g., Oral Drug Candidate")
            description = st.text_area("Description", placeholder="Description of the target product profile")

            st.markdown("**Criteria**")

            # Initialize criteria in session state
            if "tpp_criteria" not in st.session_state:
                st.session_state["tpp_criteria"] = []

            # Display existing criteria
            for i, criterion in enumerate(st.session_state["tpp_criteria"]):
                col1, col2, col3, col4, col5 = st.columns([2, 1, 1, 1, 1])
                with col1:
                    prop = st.selectbox(
                        "Property",
                        ["molecular_weight", "logp", "hbd_count", "hba_count", "rotatable_bonds"],
                        index=["molecular_weight", "logp", "hbd_count", "hba_count", "rotatable_bonds"].index(criterion.get("property", "molecular_weight")) if criterion.get("property") in ["molecular_weight", "logp", "hbd_count", "hba_count", "rotatable_bonds"] else 0,
                        key=f"prop_{i}"
                    )
                with col2:
                    min_val = st.number_input("Min", value=float(criterion.get("min", 0)), key=f"min_{i}")
                with col3:
                    max_val = st.number_input("Max", value=float(criterion.get("max", 100)), key=f"max_{i}")
                with col4:
                    weight = st.number_input("Weight", value=float(criterion.get("weight", 1.0)), min_value=0.0, step=0.1, key=f"weight_{i}")
                with col5:
                    unit = st.text_input("Unit", value=criterion.get("unit", ""), key=f"unit_{i}")

                # Update criterion
                st.session_state["tpp_criteria"][i] = {
                    "property": prop,
                    "min": min_val,
                    "max": max_val,
                    "weight": weight,
                    "unit": unit,
                }

            col1, col2 = st.columns([1, 1])
            with col1:
                if st.form_submit_button("➕ Add Criterion"):
                    st.session_state["tpp_criteria"].append({
                        "property": "molecular_weight",
                        "min": 0.0,
                        "max": 100.0,
                        "weight": 1.0,
                        "unit": "",
                    })
                    st.rerun()

            if st.session_state["tpp_criteria"]:
                with col2:
                    if st.form_submit_button("🗑️ Remove Last"):
                        if st.session_state["tpp_criteria"]:
                            st.session_state["tpp_criteria"].pop()
                            st.rerun()

            submitted = st.form_submit_button("💾 Save TPP", type="primary")

            if submitted:
                if not name:
                    st.error("TPP name is required")
                elif not st.session_state["tpp_criteria"]:
                    st.error("At least one criterion is required")
                else:
                    tpp = TargetProductProfile(
                        name=name,
                        description=description if description.strip() else None,
                        criteria=st.session_state["tpp_criteria"],
                        is_active=True,
                        created_by_id=UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None,
                    )
                    db.add(tpp)
                    db.commit()
                    st.success(f"TPP '{name}' created!")
                    st.session_state["tpp_criteria"] = []
                    st.rerun()


def render_score_tab() -> None:
    """Render the Score Compound tab."""
    st.subheader("Score Compound Against TPP")

    with db_session() as db:
        compounds = db.query(Compound).order_by(Compound.compound_id).limit(200).all()
        tpps = db.query(TargetProductProfile).filter(TargetProductProfile.is_active).all()

        if not compounds:
            st.info("No compounds available.")
            return
        if not tpps:
            st.info("No TPPs available. Create one in the 'Define TPP' tab.")
            return

        compound_options = {f"{c.compound_id} - {c.smiles[:50]}": c.id for c in compounds}
        tpp_options = {tpp.name: tpp.id for tpp in tpps}

        col1, col2 = st.columns(2)
        with col1:
            selected_compound_name = st.selectbox("Select Compound", list(compound_options.keys()))
        with col2:
            selected_tpp_name = st.selectbox("Select TPP", list(tpp_options.keys()))

        if st.button("🔍 Score", type="primary"):
            compound_id = compound_options[selected_compound_name]
            tpp_id = tpp_options[selected_tpp_name]

            with st.spinner("Scoring compound..."):
                try:
                    result = score_compound_against_tpp(compound_id, tpp_id, db)
                    st.session_state["scoring_result"] = result
                    st.session_state["scored_compound_id"] = str(compound_id)
                    st.session_state["scored_tpp_id"] = str(tpp_id)
                    st.rerun()
                except Exception as e:
                    st.error(f"Scoring failed: {e}")

        if st.session_state.get("scoring_result"):
            result = st.session_state["scoring_result"]

            # Overall Score
            overall_score = result.get("overall_score", 0)
            if overall_score >= 80:
                st.metric("Overall Score", f"{overall_score:.1f}/100", delta=None, delta_color="normal")
            elif overall_score >= 50:
                st.metric("Overall Score", f"{overall_score:.1f}/100", delta=None, delta_color="off")
            else:
                st.metric("Overall Score", f"{overall_score:.1f}/100", delta=None, delta_color="inverse")

            st.markdown("---")

            # Scorecard table
            st.markdown("### Scorecard")
            criteria_scores = result.get("criteria_scores", [])

            if criteria_scores:
                scorecard_data = []
                for criterion in criteria_scores:
                    traffic_light = criterion.get("traffic_light", "red")
                    emoji = {"green": "🟢", "yellow": "🟡", "red": "🔴"}.get(traffic_light, "⚪")

                    scorecard_data.append({
                        "Property": criterion.get("property", ""),
                        "Value": f"{criterion.get('value', 0):.2f} {criterion.get('unit', '')}",
                        "Target Range": f"{criterion.get('min', 0):.2f} - {criterion.get('max', 0):.2f}",
                        "Status": f"{emoji} {traffic_light.upper()}",
                        "Score": f"{criterion.get('score', 0):.0f}",
                    })

                df_scorecard = pd.DataFrame(scorecard_data)
                st.dataframe(df_scorecard, use_container_width=True, hide_index=True)

            # Nominate button
            st.markdown("---")
            with st.form("nominate_form"):
                notes = st.text_area("Notes (optional)", placeholder="Add notes about this nomination...")
                if st.form_submit_button("📌 Nominate Candidate", type="primary"):
                    user = get_current_user()
                    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None

                    try:
                        nomination = nominate_candidate(
                            UUID(st.session_state["scored_compound_id"]),
                            UUID(st.session_state["scored_tpp_id"]),
                            user_id,
                            notes if notes.strip() else None,
                            db
                        )
                        st.success(f"Candidate nominated! Nomination ID: {nomination.id}")
                        st.session_state.pop("scoring_result", None)
                        st.rerun()
                    except Exception as e:
                        st.error(f"Nomination failed: {e}")


def render_nominations_tab() -> None:
    """Render the Nominations tab."""
    st.subheader("Candidate Nominations")

    with db_session() as db:
        nominations = db.query(CandidateNomination).order_by(CandidateNomination.created_at.desc()).limit(100).all()

        if not nominations:
            st.info("No nominations yet.")
            return

        # Display as table
        nomination_data = []
        for nom in nominations:
            compound_name = nom.compound.compound_id if nom.compound else "Unknown"
            tpp_name = nom.tpp.name if nom.tpp else "Unknown"
            status_emoji = {"approved": "✅", "rejected": "❌", "nominated": "⏳"}.get(nom.status, "⚪")

            nomination_data.append({
                "Compound": compound_name,
                "TPP": tpp_name,
                "Score": f"{nom.overall_score:.1f}" if nom.overall_score else "N/A",
                "Status": f"{status_emoji} {nom.status}",
                "Created": nom.created_at.strftime("%Y-%m-%d") if nom.created_at else "",
            })

        df = pd.DataFrame(nomination_data)
        st.dataframe(df, use_container_width=True, hide_index=True)

        # Status update
        st.markdown("---")
        st.markdown("### Update Status")

        nom_options = {f"{nom.compound.compound_id if nom.compound else 'Unknown'} vs {nom.tpp.name if nom.tpp else 'Unknown'}": nom.id for nom in nominations}
        selected_nom_name = st.selectbox("Select Nomination", list(nom_options.keys()))

        col1, col2 = st.columns([1, 1])
        with col1:
            new_status = st.selectbox("New Status", ["nominated", "approved", "rejected"], key="nom_status_select")
        with col2:
            if st.button("Update Status", type="primary"):
                nomination = db.query(CandidateNomination).filter(CandidateNomination.id == nom_options[selected_nom_name]).first()
                if nomination:
                    nomination.status = new_status
                    db.commit()
                    st.success("Status updated!")
                    st.rerun()

        # Show details for selected nomination
        if selected_nom_name:
            nomination = db.query(CandidateNomination).filter(CandidateNomination.id == nom_options[selected_nom_name]).first()
            if nomination:
                with st.expander("Nomination Details"):
                    st.json({
                        "compound_id": str(nomination.compound_id),
                        "tpp_id": str(nomination.tpp_id),
                        "overall_score": nomination.overall_score,
                        "status": nomination.status,
                        "scores": nomination.scores,
                        "notes": nomination.notes,
                        "created_at": nomination.created_at.isoformat() if nomination.created_at else None,
                    })
