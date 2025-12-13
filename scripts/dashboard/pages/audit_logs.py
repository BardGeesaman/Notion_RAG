"""Audit Log Viewer - Admin only."""
import streamlit as st
import pandas as pd
from datetime import datetime, timedelta
from amprenta_rag.database.models import AuditLog
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.session import get_current_user


def render_audit_logs_page():
    user = get_current_user()

    # Admin only
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can view audit logs.")
        return

    st.title("📋 Audit Logs")
    st.markdown("View system activity and user actions")

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        action_filter = st.selectbox("Action", ["All", "login", "logout", "create", "update", "delete"])
    with col2:
        entity_filter = st.text_input("Entity Type", placeholder="e.g., experiment")
    with col3:
        days_back = st.selectbox("Time Range", [1, 7, 30, 90], index=1, format_func=lambda x: f"Last {x} days")

    # Query logs
    with db_session() as db:
        query = db.query(AuditLog).order_by(AuditLog.timestamp.desc())

        # Apply filters
        if action_filter != "All":
            query = query.filter(AuditLog.action == action_filter)
        if entity_filter:
            query = query.filter(AuditLog.entity_type.ilike(f"%{entity_filter}%"))

        cutoff = datetime.utcnow() - timedelta(days=days_back)
        query = query.filter(AuditLog.timestamp >= cutoff)

        logs = query.limit(500).all()

        if not logs:
            st.info("No audit logs found for the selected filters.")
            return

        # Display as table
        data = []
        for log in logs:
            data.append({
                "Timestamp": log.timestamp.strftime("%Y-%m-%d %H:%M:%S") if log.timestamp else "",
                "User": log.username or "system",
                "Action": log.action,
                "Entity": log.entity_type or "-",
                "Entity ID": (log.entity_id[:8] + "...") if log.entity_id and len(log.entity_id) > 8 else (log.entity_id or "-"),
                "Details": str(log.details)[:50] if log.details else "-"
            })

        df = pd.DataFrame(data)
        st.dataframe(df, use_container_width=True, hide_index=True)

        st.caption(f"Showing {len(logs)} most recent entries (max 500)")

        # Download option
        csv = df.to_csv(index=False)
        st.download_button("📥 Download CSV", csv, "audit_logs.csv", "text/csv")


if __name__ == "__main__":
    render_audit_logs_page()
