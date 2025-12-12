"""Generic Assay Results page for the Streamlit dashboard."""

from __future__ import annotations

import json
from typing import List, Optional
from uuid import UUID

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import GenericAssayResult, Experiment, Compound
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_generic_assays_page() -> None:
    """Render the Generic Assay Results page."""
    st.header("🧪 Generic Assay Results")
    st.markdown("Manage and view generic assay results linked to experiments or compounds.")

    tabs = st.tabs(["Browse", "Add Result"])
    
    with tabs[0]:
        _render_browse_tab()
    
    with tabs[1]:
        _render_add_tab()


def _render_browse_tab() -> None:
    """Render the browse tab with filters and results table."""
    with db_session() as db:
        # Filters
        col1, col2 = st.columns(2)
        with col1:
            # Get unique assay types
            assay_types = db.query(GenericAssayResult.assay_type).distinct().all()
            assay_type_options = ["All"] + [t[0] for t in assay_types if t[0]]
            selected_assay_type = st.selectbox("Filter by Assay Type", assay_type_options)
        
        with col2:
            # Get experiments that have assay results
            # Get distinct experiment IDs first, then fetch experiments
            # Avoid DISTINCT on Experiment rows which have JSON columns
            exp_ids = db.query(GenericAssayResult.experiment_id).filter(
                GenericAssayResult.experiment_id.isnot(None)
            ).distinct().all()
            exp_ids = [e[0] for e in exp_ids if e[0]]
            experiments_with_results = db.query(Experiment).filter(Experiment.id.in_(exp_ids)).all() if exp_ids else []
            experiment_options = ["All"] + [exp.name for exp in experiments_with_results]
            selected_experiment = st.selectbox("Filter by Experiment", experiment_options)
        
        # Build query
        query = db.query(GenericAssayResult)
        
        if selected_assay_type != "All":
            query = query.filter(GenericAssayResult.assay_type == selected_assay_type)
        
        if selected_experiment != "All":
            exp = db.query(Experiment).filter(Experiment.name == selected_experiment).first()
            if exp:
                query = query.filter(GenericAssayResult.experiment_id == exp.id)
        
        results = query.order_by(GenericAssayResult.created_at.desc()).all()
        
        st.metric("Total Results", len(results))
        
        if results:
            # Display results table
            result_data = []
            for result in results:
                experiment_name = ""
                compound_id = ""
                if result.experiment_id:
                    exp = db.query(Experiment).filter(Experiment.id == result.experiment_id).first()
                    experiment_name = exp.name if exp else ""
                if result.compound_id:
                    compound_id = str(result.compound_id)[:8] + "..."
                
                result_data.append({
                    "Assay Name": result.assay_name,
                    "Assay Type": result.assay_type,
                    "Experiment": experiment_name or "N/A",
                    "Compound ID": compound_id or "N/A",
                    "Created": result.created_at.strftime("%Y-%m-%d %H:%M") if result.created_at else "",
                })
            
            df = pd.DataFrame(result_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Detail view
            st.markdown("---")
            st.subheader("Result Details")
            
            selected_result_name = st.selectbox(
                "Select result to view details",
                [r.assay_name for r in results],
                key="result_detail_select"
            )
            
            if selected_result_name:
                result = next(r for r in results if r.assay_name == selected_result_name)
                
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**Assay Name:** {result.assay_name}")
                    st.write(f"**Assay Type:** {result.assay_type}")
                    if result.experiment_id:
                        exp = db.query(Experiment).filter(Experiment.id == result.experiment_id).first()
                        st.write(f"**Experiment:** {exp.name if exp else 'N/A'}")
                    if result.compound_id:
                        st.write(f"**Compound ID:** `{result.compound_id}`")
                    st.write(f"**Created:** {result.created_at.strftime('%Y-%m-%d %H:%M') if result.created_at else 'N/A'}")
                
                with col2:
                    if result.created_by_id:
                        from amprenta_rag.database.models import User
                        user = db.query(User).filter(User.id == result.created_by_id).first()
                        st.write(f"**Created By:** {user.username if user else 'Unknown'}")
                
                st.markdown("---")
                st.subheader("Result Data")
                st.json(result.result_data)
                
                if result.assay_metadata:
                    st.subheader("Metadata")
                    st.json(result.assay_metadata)
        else:
            st.info("No results found. Add a new result using the 'Add Result' tab.")


def _render_add_tab() -> None:
    """Render the add result tab with form."""
    user = get_current_user()
    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None
    
    if not user_id:
        st.warning("Please log in to add assay results.")
        return
    
    with db_session() as db:
        st.subheader("Add New Assay Result")
        
        with st.form("add_assay_result", clear_on_submit=True):
            assay_name = st.text_input("Assay Name*", placeholder="e.g., EGFR IC50 Assay")
            assay_type = st.selectbox(
                "Assay Type*",
                ["biochemical", "cellular", "in_vivo", "custom", "other"]
            )
            
            # Link to experiment or compound
            link_type = st.radio("Link to", ["Experiment", "Compound", "None"], horizontal=True)
            
            experiment_id = None
            compound_id = None
            
            if link_type == "Experiment":
                experiments = db.query(Experiment).order_by(Experiment.name).all()
                if experiments:
                    exp_options = {exp.name: exp.id for exp in experiments}
                    selected_exp_name = st.selectbox("Select Experiment", list(exp_options.keys()))
                    experiment_id = exp_options[selected_exp_name]
                else:
                    st.info("No experiments available.")
            
            elif link_type == "Compound":
                compounds = db.query(Compound).order_by(Compound.compound_id).limit(100).all()
                if compounds:
                    compound_options = {f"{c.compound_id} ({c.smiles[:30]}...)" if len(c.smiles) > 30 else f"{c.compound_id} ({c.smiles})": c.id for c in compounds}
                    selected_compound_label = st.selectbox("Select Compound", list(compound_options.keys()))
                    compound_id = compound_options[selected_compound_label]
                else:
                    st.info("No compounds available.")
            
            # JSON editor for result_data
            st.markdown("**Result Data (JSON)*:**")
            st.caption("Enter the assay result data as JSON. Example: {\"value\": 10.5, \"unit\": \"nM\", \"replicate_count\": 3}")
            result_data_text = st.text_area(
                "Result Data",
                placeholder='{"value": 10.5, "unit": "nM"}',
                height=150,
                key="result_data_input"
            )
            
            # Optional metadata
            st.markdown("**Metadata (JSON, optional):**")
            st.caption("Additional metadata such as conditions, notes, etc.")
            metadata_text = st.text_area(
                "Metadata",
                placeholder='{"conditions": {"temperature": "37C", "ph": 7.4}, "notes": "..."}',
                height=100,
                key="metadata_input"
            )
            
            submitted = st.form_submit_button("💾 Save Result", type="primary")
            
            if submitted:
                if not assay_name or not assay_name.strip():
                    st.error("Assay name is required.")
                elif not result_data_text or not result_data_text.strip():
                    st.error("Result data is required.")
                else:
                    try:
                        # Parse JSON
                        result_data = json.loads(result_data_text)
                        assay_metadata = json.loads(metadata_text) if metadata_text.strip() else None
                        
                        # Create result
                        result = GenericAssayResult(
                            assay_name=assay_name.strip(),
                            assay_type=assay_type,
                            experiment_id=experiment_id,
                            compound_id=compound_id,
                            result_data=result_data,
                            assay_metadata=assay_metadata,
                            created_by_id=user_id,
                        )
                        
                        db.add(result)
                        db.commit()
                        st.success(f"Assay result '{assay_name}' added successfully!")
                        st.rerun()
                        
                    except json.JSONDecodeError as e:
                        st.error(f"Invalid JSON format: {e}")
                    except Exception as e:
                        st.error(f"Failed to add assay result: {e}")


__all__ = ["render_generic_assays_page"]

