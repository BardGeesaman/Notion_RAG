"""Feedback and feature requests page."""
import streamlit as st
from datetime import datetime
from uuid import UUID
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Feedback
from amprenta_rag.auth.session import get_current_user


def render_feedback_page():
    st.title("💬 Feedback & Feature Requests")
    st.markdown("Submit feedback, report bugs, or request features")

    tab1, tab2, tab3 = st.tabs(["Submit Feedback", "My Submissions", "Admin Triage"])

    with tab1:
        render_submit_feedback_tab()

    with tab2:
        render_my_submissions_tab()

    with tab3:
        render_admin_triage_tab()


def render_submit_feedback_tab():
    st.subheader("Submit Feedback")
    user = get_current_user()

    with st.form("submit_feedback"):
        feedback_type = st.selectbox("Type*", ["Bug", "Feature"])
        title = st.text_input("Title*", max_chars=255)
        description = st.text_area("Description", height=150, placeholder="Describe the bug or feature request in detail...")

        submitted = st.form_submit_button("Submit Feedback", type="primary")

        if submitted:
            if not title:
                st.error("Title is required")
                return

            db_gen = get_db()
            db = next(db_gen)
            try:
                feedback = Feedback(
                    user_id=UUID(user.get("id")) if user and user.get("id") != "test" else None,
                    feedback_type=feedback_type.lower(),
                    title=title,
                    description=description if description else None,
                    status="new",
                    priority="medium",
                )
                db.add(feedback)
                db.commit()
                st.success(f"Feedback submitted! Thank you for your input.")
                st.rerun()
            finally:
                db_gen.close()


def render_my_submissions_tab():
    st.subheader("My Submissions")
    user = get_current_user()

    if not user:
        st.error("Please log in to view your submissions")
        return

    db_gen = get_db()
    db = next(db_gen)
    try:
        user_id = user.get("id")
        if user_id == "test":
            st.info("Test user - no submissions to display")
            return

        feedbacks = db.query(Feedback).filter(
            Feedback.user_id == UUID(user_id)
        ).order_by(Feedback.created_at.desc()).all()

        if not feedbacks:
            st.info("You haven't submitted any feedback yet.")
            return

        for fb in feedbacks:
            status_colors = {
                "new": "🆕",
                "reviewing": "👀",
                "planned": "📋",
                "done": "✅",
                "wontfix": "❌",
            }
            status_icon = status_colors.get(fb.status, "❓")

            priority_colors = {
                "low": "🟢",
                "medium": "🟡",
                "high": "🟠",
                "critical": "🔴",
            }
            priority_icon = priority_colors.get(fb.priority, "⚪")

            with st.expander(f"{status_icon} {priority_icon} **{fb.title}** ({fb.feedback_type})"):
                st.markdown(f"**Status:** {fb.status}")
                st.markdown(f"**Priority:** {fb.priority}")
                if fb.description:
                    st.markdown(f"**Description:** {fb.description}")
                st.caption(f"Submitted: {fb.created_at.strftime('%Y-%m-%d %H:%M') if fb.created_at else 'Unknown'}")
                if fb.updated_at and fb.updated_at != fb.created_at:
                    st.caption(f"Updated: {fb.updated_at.strftime('%Y-%m-%d %H:%M')}")
    finally:
        db_gen.close()


def render_admin_triage_tab():
    st.subheader("Admin Triage")
    user = get_current_user()

    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can view this page.")
        return

    db_gen = get_db()
    db = next(db_gen)
    try:
        feedbacks = db.query(Feedback).order_by(Feedback.created_at.desc()).all()

        if not feedbacks:
            st.info("No feedback submissions yet.")
            return

        st.metric("Total Submissions", len(feedbacks))

        for fb in feedbacks:
            status_colors = {
                "new": "🆕",
                "reviewing": "👀",
                "planned": "📋",
                "done": "✅",
                "wontfix": "❌",
            }
            status_icon = status_colors.get(fb.status, "❓")

            priority_colors = {
                "low": "🟢",
                "medium": "🟡",
                "high": "🟠",
                "critical": "🔴",
            }
            priority_icon = priority_colors.get(fb.priority, "⚪")

            with st.expander(f"{status_icon} {priority_icon} **{fb.title}** ({fb.feedback_type})"):
                col1, col2 = st.columns(2)
                with col1:
                    new_status = st.selectbox(
                        "Status",
                        ["new", "reviewing", "planned", "done", "wontfix"],
                        index=["new", "reviewing", "planned", "done", "wontfix"].index(fb.status) if fb.status in ["new", "reviewing", "planned", "done", "wontfix"] else 0,
                        key=f"status_{fb.id}",
                    )
                with col2:
                    new_priority = st.selectbox(
                        "Priority",
                        ["low", "medium", "high", "critical"],
                        index=["low", "medium", "high", "critical"].index(fb.priority) if fb.priority in ["low", "medium", "high", "critical"] else 1,
                        key=f"priority_{fb.id}",
                    )

                if st.button("Update", key=f"update_{fb.id}"):
                    fb.status = new_status
                    fb.priority = new_priority
                    fb.updated_at = datetime.utcnow()
                    db.commit()
                    st.success("Updated!")
                    st.rerun()

                st.markdown("---")
                if fb.description:
                    st.markdown(f"**Description:** {fb.description}")
                st.caption(f"Submitted: {fb.created_at.strftime('%Y-%m-%d %H:%M') if fb.created_at else 'Unknown'}")
                if fb.user_id:
                    # Could fetch username if needed
                    st.caption(f"User ID: {str(fb.user_id)[:8]}...")
    finally:
        db_gen.close()


if __name__ == "__main__":
    render_feedback_page()
