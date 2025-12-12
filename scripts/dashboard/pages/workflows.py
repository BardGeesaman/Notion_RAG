"""Workflow Automation Management page."""
from __future__ import annotations

import json
from datetime import datetime

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import WorkflowRule, WorkflowExecution
from amprenta_rag.automation.engine import TRIGGER_TYPES
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_workflows_page() -> None:
    """Render the Workflows page."""
    st.header("⚙️ Workflow Automation")
    st.markdown("Create and manage automated workflow rules.")
    
    tab1, tab2, tab3 = st.tabs(["Rules", "Create Rule", "History"])
    
    with tab1:
        render_rules_tab()
    with tab2:
        render_create_tab()
    with tab3:
        render_history_tab()


def render_rules_tab() -> None:
    """Render the Rules tab."""
    st.subheader("Workflow Rules")
    
    with db_session() as db:
        rules = db.query(WorkflowRule).order_by(WorkflowRule.created_at.desc()).all()
        
        if not rules:
            st.info("No workflow rules yet. Create your first rule!")
            return
        
        for rule in rules:
            with st.expander(f"**{rule.name}** - {rule.trigger_type} → {rule.action_type}"):
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.markdown(f"**Description:** {rule.description or 'No description'}")
                    st.markdown(f"**Status:** {'✅ Active' if rule.is_active else '❌ Inactive'}")
                    if rule.trigger_config:
                        st.markdown(f"**Trigger Config:**")
                        st.json(rule.trigger_config)
                    if rule.action_config:
                        st.markdown(f"**Action Config:**")
                        st.json(rule.action_config)
                with col2:
                    # Toggle active status
                    new_status = st.checkbox(
                        "Active",
                        value=rule.is_active,
                        key=f"rule_active_{rule.id}"
                    )
                    if new_status != rule.is_active:
                        rule.is_active = new_status
                        db.commit()
                        st.success("Rule status updated!")
                        st.rerun()
                    
                    if st.button("🗑️ Delete", key=f"delete_rule_{rule.id}"):
                        db.delete(rule)
                        db.commit()
                        st.success("Rule deleted!")
                        st.rerun()


def render_create_tab() -> None:
    """Render the Create Rule tab."""
    st.subheader("Create New Workflow Rule")
    
    user = get_current_user()
    if not user:
        st.error("You must be logged in to create workflow rules.")
        return
    
    with db_session() as db:
        with st.form("create_workflow_rule"):
            name = st.text_input("Rule Name*", placeholder="e.g., Notify on new experiment")
            description = st.text_area("Description", placeholder="Optional description")
            
            trigger_type = st.selectbox("Trigger Type*", TRIGGER_TYPES)
            
            action_type = st.selectbox(
                "Action Type*",
                ["send_notification", "add_note", "run_validation"]
            )
            
            # Action config as JSON
            st.markdown("**Action Configuration (JSON)**")
            action_config_json = st.text_area(
                "Action Config",
                placeholder='{"user_id": "creator", "title": "New {entity_type} created"}',
                height=150,
                help="JSON object with action parameters"
            )
            
            # Trigger config (optional)
            st.markdown("**Trigger Configuration (Optional JSON)**")
            trigger_config_json = st.text_area(
                "Trigger Config",
                placeholder='{"design_type": "case_control"}',
                height=100,
                help="JSON object with filter conditions - leave empty to match all"
            )
            
            submitted = st.form_submit_button("Create Rule", type="primary")
            
            if submitted:
                if not name or not trigger_type or not action_type:
                    st.error("Name, trigger type, and action type are required")
                    return
                
                # Parse action config
                try:
                    action_config = json.loads(action_config_json) if action_config_json.strip() else {}
                except json.JSONDecodeError as e:
                    st.error(f"Invalid action config JSON: {e}")
                    return
                
                # Parse trigger config
                trigger_config = None
                if trigger_config_json.strip():
                    try:
                        trigger_config = json.loads(trigger_config_json)
                    except json.JSONDecodeError as e:
                        st.error(f"Invalid trigger config JSON: {e}")
                        return
                
                # Create rule
                rule = WorkflowRule(
                    name=name,
                    description=description if description.strip() else None,
                    trigger_type=trigger_type,
                    trigger_config=trigger_config,
                    action_type=action_type,
                    action_config=action_config,
                    is_active=True,
                    created_by_id=user.get("id") if user.get("id") and user.get("id") != "test" else None,
                )
                
                db.add(rule)
                db.commit()
                st.success(f"Workflow rule '{name}' created!")
                st.rerun()


def render_history_tab() -> None:
    """Render the History tab."""
    st.subheader("Execution History")
    
    with db_session() as db:
        executions = db.query(WorkflowExecution).order_by(
            WorkflowExecution.triggered_at.desc()
        ).limit(100).all()
        
        if not executions:
            st.info("No workflow executions yet.")
            return
        
        # Display as table
        execution_data = []
        for exec in executions:
            rule_name = exec.rule.name if exec.rule else "Unknown"
            status_emoji = "✅" if exec.status == "success" else "❌" if exec.status == "failed" else "⏳"
            triggered_at_str = exec.triggered_at.strftime("%Y-%m-%d %H:%M:%S") if exec.triggered_at else "Unknown"
            
            execution_data.append({
                "Rule": rule_name,
                "Status": f"{status_emoji} {exec.status}",
                "Triggered": triggered_at_str,
            })
        
        df = pd.DataFrame(execution_data)
        st.dataframe(df, use_container_width=True, hide_index=True)
