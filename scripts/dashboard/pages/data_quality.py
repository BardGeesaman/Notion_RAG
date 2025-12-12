"""Data Quality page for validation and quality checks."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.utils.validation import run_all_validations, ValidationIssue
from scripts.dashboard.db_session import db_session


def render_data_quality_page() -> None:
    """Render the Data Quality page."""
    st.header("✅ Data Quality")
    st.markdown("Run validation checks on experiments and compounds to identify data quality issues.")
    
    if st.button("🔍 Run Validation", type="primary"):
        with st.spinner("Running validation checks..."):
            with db_session() as db:
                issues = run_all_validations(db)
        
        # Store in session state
        st.session_state["validation_issues"] = issues
    
    # Display results if available
    if "validation_issues" in st.session_state:
        issues = st.session_state["validation_issues"]
        
        # Summary stats
        error_count = sum(1 for i in issues if i.severity == "error")
        warning_count = sum(1 for i in issues if i.severity == "warning")
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Errors", error_count, delta=None, delta_color="inverse")
        with col2:
            st.metric("Warnings", warning_count, delta=None)
        
        st.markdown("---")
        
        if issues:
            # Convert to DataFrame for display
            issues_data = []
            for issue in issues:
                issues_data.append({
                    "Entity Type": issue.entity_type,
                    "Entity ID": issue.entity_id[:8] + "..." if len(issue.entity_id) > 8 else issue.entity_id,
                    "Field": issue.field,
                    "Issue": issue.issue,
                    "Severity": issue.severity,
                })
            
            df = pd.DataFrame(issues_data)
            
            # Color code by severity
            def color_severity(val):
                if val == "error":
                    return "background-color: #ffcccc"
                elif val == "warning":
                    return "background-color: #fff4cc"
                return ""
            
            styled_df = df.style.map(color_severity, subset=["Severity"])
            st.dataframe(styled_df, use_container_width=True, hide_index=True)
            
            # Download option
            csv = df.to_csv(index=False)
            st.download_button(
                label="📥 Download Issues (CSV)",
                data=csv,
                file_name="validation_issues.csv",
                mime="text/csv",
            )
        else:
            st.success("✅ No validation issues found! All data passes validation checks.")
    else:
        st.info("👆 Click 'Run Validation' to check data quality.")
