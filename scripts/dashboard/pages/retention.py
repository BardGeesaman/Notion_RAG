"""Data Retention page for managing retention policies and archived entities."""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional
from uuid import UUID

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import RetentionPolicy, Experiment, Dataset
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_retention_page() -> None:
    """Render the Data Retention page."""
    user = get_current_user()
    
    # Admin only
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can manage data retention.")
        return
    
    st.header("🗄️ Data Retention")
    st.markdown("Manage retention policies and archived entities.")

    tabs = st.tabs(["Policies", "Archived"])
    
    with tabs[0]:
        _render_policies_tab()
    
    with tabs[1]:
        _render_archived_tab()


def _render_policies_tab() -> None:
    """Render the policies tab."""
    with db_session() as db:
        st.subheader("Retention Policies")
        
        # Get all policies
        policies = db.query(RetentionPolicy).order_by(RetentionPolicy.created_at.desc()).all()
        
        if policies:
            st.metric("Total Policies", len(policies))
            
            # Display policies table
            policy_data = []
            for policy in policies:
                policy_data.append({
                    "Name": policy.name,
                    "Entity Type": policy.entity_type,
                    "Retention Days": policy.retention_days,
                    "Action": policy.action,
                    "Status": "✅ Active" if policy.is_active else "❌ Inactive",
                    "Created": policy.created_at.strftime("%Y-%m-%d") if policy.created_at else "",
                })
            
            df_policies = pd.DataFrame(policy_data)
            st.dataframe(df_policies, use_container_width=True, hide_index=True)
            
            # Toggle active status
            st.markdown("---")
            st.subheader("Toggle Policy Status")
            policy_options = {f"{p.name} ({p.entity_type})": p.id for p in policies}
            if policy_options:
                selected_policy_label = st.selectbox("Select Policy", list(policy_options.keys()))
                selected_policy_id = policy_options[selected_policy_label]
                
                policy = db.query(RetentionPolicy).filter(RetentionPolicy.id == selected_policy_id).first()
                if policy:
                    current_status = "Active" if policy.is_active else "Inactive"
                    new_status = "Inactive" if policy.is_active else "Active"
                    
                    if st.button(f"Toggle to {new_status}", key=f"toggle_{selected_policy_id}"):
                        policy.is_active = not policy.is_active
                        db.commit()
                        st.success(f"Policy '{policy.name}' is now {new_status.lower()}.")
                        st.rerun()
        else:
            st.info("No retention policies defined yet.")
        
        st.markdown("---")
        st.subheader("Create New Policy")
        
        with st.form("create_policy", clear_on_submit=True):
            policy_name = st.text_input("Policy Name*", placeholder="e.g., Archive Old Experiments")
            entity_type = st.selectbox(
                "Entity Type*",
                ["experiment", "dataset", "compound", "signature", "program"]
            )
            retention_days = st.number_input("Retention Days*", min_value=1, value=365, step=1)
            action = st.selectbox(
                "Action*",
                ["archive", "delete", "notify"]
            )
            is_active = st.checkbox("Active", value=True)
            
            submitted = st.form_submit_button("💾 Create Policy", type="primary")
            
            if submitted:
                if not policy_name or not policy_name.strip():
                    st.error("Policy name is required.")
                else:
                    try:
                        policy = RetentionPolicy(
                            name=policy_name.strip(),
                            entity_type=entity_type,
                            retention_days=retention_days,
                            action=action,
                            is_active=is_active,
                        )
                        
                        db.add(policy)
                        db.commit()
                        st.success(f"Policy '{policy_name}' created successfully!")
                        st.rerun()
                        
                    except Exception as e:
                        st.error(f"Failed to create policy: {e}")


def _render_archived_tab() -> None:
    """Render the archived entities tab."""
    with db_session() as db:
        st.subheader("Archived Entities")
        
        # Entity type filter
        entity_type_filter = st.radio(
            "Entity Type",
            ["All", "Experiments", "Datasets"],
            horizontal=True
        )
        
        # Get archived entities
        archived_experiments = []
        archived_datasets = []
        
        if entity_type_filter in ["All", "Experiments"]:
            archived_experiments = db.query(Experiment).filter(Experiment.is_archived == True).order_by(Experiment.archived_at.desc()).all()
        
        if entity_type_filter in ["All", "Datasets"]:
            archived_datasets = db.query(Dataset).filter(Dataset.is_archived == True).order_by(Dataset.archived_at.desc()).all()
        
        total_archived = len(archived_experiments) + len(archived_datasets)
        st.metric("Total Archived", total_archived)
        
        if archived_experiments or archived_datasets:
            # Display archived experiments
            if archived_experiments:
                st.markdown("### Archived Experiments")
                exp_data = []
                for exp in archived_experiments:
                    exp_data.append({
                        "Name": exp.name,
                        "Type": exp.type or "N/A",
                        "Archived At": exp.archived_at.strftime("%Y-%m-%d %H:%M") if exp.archived_at else "N/A",
                        "Created": exp.created_at.strftime("%Y-%m-%d") if exp.created_at else "",
                    })
                
                df_exps = pd.DataFrame(exp_data)
                
                # Add restore buttons
                for idx, exp in enumerate(archived_experiments):
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        st.write(f"**{exp.name}**")
                        st.caption(f"Archived: {exp.archived_at.strftime('%Y-%m-%d %H:%M') if exp.archived_at else 'N/A'}")
                    with col2:
                        if st.button("Restore", key=f"restore_exp_{exp.id}"):
                            exp.is_archived = False
                            exp.archived_at = None
                            db.commit()
                            st.success(f"Experiment '{exp.name}' restored.")
                            st.rerun()
                    st.divider()
            
            # Display archived datasets
            if archived_datasets:
                st.markdown("### Archived Datasets")
                for dataset in archived_datasets:
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        st.write(f"**{dataset.name}**")
                        st.caption(f"Archived: {dataset.archived_at.strftime('%Y-%m-%d %H:%M') if dataset.archived_at else 'N/A'}")
                    with col2:
                        if st.button("Restore", key=f"restore_dataset_{dataset.id}"):
                            dataset.is_archived = False
                            dataset.archived_at = None
                            db.commit()
                            st.success(f"Dataset '{dataset.name}' restored.")
                            st.rerun()
                    st.divider()
        else:
            st.info("No archived entities found.")
        
        st.markdown("---")
        st.subheader("Archive Entities")
        
        archive_entity_type = st.selectbox("Entity Type", ["Experiment", "Dataset"], key="archive_entity_type")
        
        if archive_entity_type == "Experiment":
            experiments = db.query(Experiment).filter(Experiment.is_archived == False).order_by(Experiment.name).limit(100).all()
            if experiments:
                exp_options = {exp.name: exp.id for exp in experiments}
                selected_exp_name = st.selectbox("Select Experiment", list(exp_options.keys()))
                
                if st.button("Archive Selected Experiment", type="primary"):
                    exp_id = exp_options[selected_exp_name]
                    exp = db.query(Experiment).filter(Experiment.id == exp_id).first()
                    if exp:
                        exp.is_archived = True
                        exp.archived_at = datetime.utcnow()
                        db.commit()
                        st.success(f"Experiment '{exp.name}' archived.")
                        st.rerun()
            else:
                st.info("No unarchived experiments available.")
        
        elif archive_entity_type == "Dataset":
            datasets = db.query(Dataset).filter(Dataset.is_archived == False).order_by(Dataset.name).limit(100).all()
            if datasets:
                dataset_options = {ds.name: ds.id for ds in datasets}
                selected_dataset_name = st.selectbox("Select Dataset", list(dataset_options.keys()))
                
                if st.button("Archive Selected Dataset", type="primary"):
                    dataset_id = dataset_options[selected_dataset_name]
                    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                    if dataset:
                        dataset.is_archived = True
                        dataset.archived_at = datetime.utcnow()
                        db.commit()
                        st.success(f"Dataset '{dataset.name}' archived.")
                        st.rerun()
            else:
                st.info("No unarchived datasets available.")


__all__ = ["render_retention_page"]

