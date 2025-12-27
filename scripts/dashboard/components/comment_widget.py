"""Comment widget component for displaying and adding contextual comments."""

from __future__ import annotations

import re
from uuid import UUID
import streamlit as st
from amprenta_rag.utils.comments import add_comment, get_comments, delete_comment, update_comment
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_comments_widget(entity_type: str, entity_id: str | UUID) -> None:
    """
    Render a comments widget for an entity.

    Args:
        entity_type: Type of entity ("experiment", "dataset", "signature", "compound")
        entity_id: UUID string or UUID object of the entity
    """
    # Convert entity_id to UUID if it's a string
    if isinstance(entity_id, str):
        entity_id = UUID(entity_id)

    user = get_current_user()
    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None

    with db_session() as db:
        # Get existing comments
        comments = get_comments(entity_type, entity_id, db)

        st.markdown("### 💬 Comments")

        # Display existing comments
        if comments:
            for comment in comments:
                _render_comment(comment, user_id, db, entity_type, entity_id)
        else:
            st.info("No comments yet. Be the first to comment!")

        st.markdown("---")

        # Add new comment form
        if user_id:
            with st.form(f"add_comment_{entity_type}_{entity_id}", clear_on_submit=True):
                comment_text = st.text_area("Add a comment", placeholder="Write your comment here...", height=100)
                submitted = st.form_submit_button("💬 Add Comment", type="primary")

                if submitted:
                    if comment_text.strip():
                        try:
                            add_comment(
                                entity_type=entity_type,
                                entity_id=entity_id,
                                content=comment_text.strip(),
                                user_id=user_id,
                                db=db,
                            )
                            st.success("Comment added!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to add comment: {e}")
                    else:
                        st.warning("Please enter a comment.")
        else:
            st.info("Please log in to add comments.")


def _highlight_mentions(content: str) -> str:
    """Highlight @mentions in comment content with styled HTML.
    
    Args:
        content: Comment text
        
    Returns:
        HTML string with highlighted mentions
    """
    pattern = r'@([a-zA-Z0-9_-]+)'
    highlighted = re.sub(
        pattern,
        r'<span style="color: #1f77b4; font-weight: bold; background-color: #e3f2fd; padding: 2px 4px; border-radius: 3px;">@\1</span>',
        content
    )
    return highlighted


def _render_comment(comment: dict, user_id: UUID | None, db, entity_type: str, entity_id: UUID) -> None:
    """Render a single comment with its replies."""
    comment_id = UUID(comment["id"])
    is_owner = user_id and comment.get("author_id") == str(user_id)
    is_editing = st.session_state.get(f"edit_{comment_id}", False)

    # Comment container
    with st.container():
        # Show edit form if editing
        if is_editing:
            with st.form(f"edit_form_{comment_id}"):
                edited_text = st.text_area(
                    "Edit comment",
                    value=comment["content"],
                    key=f"edit_text_{comment_id}"
                )
                col1, col2 = st.columns([1, 1])
                with col1:
                    if st.form_submit_button("💾 Save", type="primary"):
                        if edited_text.strip():
                            try:
                                updated = update_comment(comment_id, edited_text.strip(), user_id, db)
                                if updated:
                                    st.session_state[f"edit_{comment_id}"] = False
                                    st.success("Comment updated!")
                                    st.rerun()
                                else:
                                    st.error("Failed to update comment.")
                            except Exception as e:
                                st.error(f"Error: {e}")
                        else:
                            st.warning("Comment cannot be empty.")
                with col2:
                    if st.form_submit_button("Cancel"):
                        st.session_state[f"edit_{comment_id}"] = False
                        st.rerun()
        else:
            col1, col2 = st.columns([5, 1])

            with col1:
                # Show edited badge if comment was updated
                edited_badge = ""
                if comment.get("updated_at") and comment.get("created_at"):
                    if comment["updated_at"] != comment["created_at"]:
                        edited_badge = " <span style='color: gray; font-size: 0.9em;'>(edited)</span>"
                
                st.markdown(
                    f"**{comment['author']}** - {comment['created_at'].strftime('%Y-%m-%d %H:%M')}{edited_badge}",
                    unsafe_allow_html=True
                )
                
                # Highlight @mentions in content
                highlighted_content = _highlight_mentions(comment["content"])
                st.markdown(highlighted_content, unsafe_allow_html=True)

            with col2:
                # Edit and delete buttons (only for own comments)
                if is_owner:
                    if st.button("✏️", key=f"edit_{comment_id}_btn", help="Edit comment"):
                        st.session_state[f"edit_{comment_id}"] = True
                        st.rerun()
                    
                    if st.button("🗑️", key=f"delete_{comment_id}", help="Delete comment"):
                        try:
                            if delete_comment(comment_id, user_id, db):
                                st.success("Comment deleted!")
                                st.rerun()
                            else:
                                st.error("Failed to delete comment.")
                        except Exception as e:
                            st.error(f"Error: {e}")

        # Replies section
        if comment.get("replies"):
            with st.container():
                st.markdown("---")
                for reply in comment["replies"]:
                    st.markdown(f"↳ **{reply.get('author', 'Unknown')}** - {reply['created_at'].strftime('%Y-%m-%d %H:%M')}")
                    st.markdown(f"  {reply['content']}")

        # Reply form
        if user_id:
            if st.session_state.get(f"show_reply_{comment_id}", False):
                with st.form(f"reply_form_{comment_id}", clear_on_submit=True):
                    reply_text = st.text_area("Write a reply", placeholder="Your reply...", key=f"reply_text_{comment_id}")
                    col1, col2 = st.columns([1, 1])
                    with col1:
                        if st.form_submit_button("💬 Reply", type="primary"):
                            if reply_text.strip():
                                try:
                                    add_comment(
                                        entity_type=entity_type,
                                        entity_id=entity_id,
                                        content=reply_text.strip(),
                                        user_id=user_id,
                                        db=db,
                                        parent_id=comment_id,
                                    )
                                    st.session_state[f"show_reply_{comment_id}"] = False
                                    st.success("Reply added!")
                                    st.rerun()
                                except Exception as e:
                                    st.error(f"Failed to add reply: {e}")
                            else:
                                st.warning("Please enter a reply.")
                    with col2:
                        if st.form_submit_button("Cancel"):
                            st.session_state[f"show_reply_{comment_id}"] = False
                            st.rerun()
            else:
                if st.button("↳ Reply", key=f"show_reply_{comment_id}"):
                    st.session_state[f"show_reply_{comment_id}"] = True
                    st.rerun()

        st.markdown("---")

