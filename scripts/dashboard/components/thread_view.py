"""Thread discussion view component for notebook reviews."""

from datetime import datetime
from typing import Dict, List, Optional
from uuid import uuid4

import httpx
import streamlit as st


def render_thread_view(review_id: str, api_base: str = "http://localhost:8000") -> None:
    """Render the thread discussion view for a notebook review.
    
    Args:
        review_id: UUID of the notebook review
        api_base: Base URL for the API server
    """
    st.subheader("💬 Discussion Threads")
    
    # Load existing threads
    try:
        with st.spinner("Loading discussion threads..."):
            response = httpx.get(f"{api_base}/api/v1/reviews/{review_id}/threads")
            response.raise_for_status()
            threads = response.json()
    except httpx.ConnectError:
        st.warning("⚠️ API server unavailable. Showing demo data.")
        threads = _get_demo_threads()
    except Exception as e:
        st.error(f"Failed to load threads: {e}")
        threads = []
    
    # New thread creation section
    st.markdown("### 📝 Start New Discussion")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        new_thread_title = st.text_input(
            "Thread Title",
            placeholder="e.g., 'Question about methodology in cell 3'",
            key="new_thread_title"
        )
    with col2:
        cell_index = st.number_input(
            "Cell # (optional)",
            min_value=0,
            value=None,
            step=1,
            key="new_thread_cell"
        )
    
    if st.button("💬 Create Thread", key="create_thread", type="primary"):
        if new_thread_title.strip():
            try:
                with st.spinner("Creating thread..."):
                    payload = {"title": new_thread_title.strip()}
                    if cell_index is not None:
                        payload["cell_index"] = int(cell_index)
                    
                    response = httpx.post(
                        f"{api_base}/api/v1/reviews/{review_id}/threads",
                        json=payload
                    )
                    response.raise_for_status()
                    
                st.success("✅ Thread created successfully!")
                st.rerun()
            except httpx.ConnectError:
                st.warning("⚠️ API server unavailable. Thread creation disabled.")
            except Exception as e:
                st.error(f"Failed to create thread: {e}")
        else:
            st.error("Please enter a thread title.")
    
    st.divider()
    
    # Display existing threads
    if not threads:
        st.info("💭 No discussions yet. Start the first thread above!")
        return
    
    st.markdown(f"### 📋 Threads ({len(threads)})")
    
    for thread in threads:
        _render_thread(thread, api_base)


def _render_thread(thread: Dict, api_base: str) -> None:
    """Render an individual thread with comments and actions."""
    thread_id = thread["id"]
    
    # Status badge
    status = thread.get("status", "open")
    status_emoji = {"open": "🔓", "resolved": "✅", "wontfix": "❌"}.get(status, "❓")
    status_color = {"open": "blue", "resolved": "green", "wontfix": "red"}.get(status, "gray")
    
    # Cell anchor display
    cell_display = ""
    if thread.get("cell_index") is not None:
        cell_display = f" • Cell #{thread['cell_index']}"
    
    # Thread header
    with st.expander(
        f"{status_emoji} **{thread['title']}**{cell_display}",
        expanded=(status == "open")
    ):
        # Thread metadata
        created_at = thread.get("created_at", "")
        if created_at:
            try:
                dt = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
                created_str = dt.strftime("%Y-%m-%d %H:%M")
            except (ValueError, AttributeError):
                created_str = created_at
        else:
            created_str = "Unknown"
        
        st.caption(f"Created: {created_str} | Status: :{status_color}[{status.upper()}]")
        
        # Thread actions
        col1, col2, col3 = st.columns([1, 1, 3])
        
        with col1:
            if status == "open":
                if st.button("✅ Resolve", key=f"resolve_{thread_id}"):
                    _update_thread_status(thread_id, "resolved", api_base)
            else:
                if st.button("🔓 Reopen", key=f"reopen_{thread_id}"):
                    _update_thread_status(thread_id, "open", api_base)
        
        with col2:
            if status == "open":
                if st.button("❌ Won't Fix", key=f"wontfix_{thread_id}"):
                    _update_thread_status(thread_id, "wontfix", api_base)
        
        # Comments section
        comments = thread.get("comments", [])
        if comments:
            st.markdown("**Comments:**")
            for comment in comments:
                _render_comment(comment, thread_id, api_base)
        else:
            st.info("💬 No comments yet. Be the first to comment!")
        
        # New comment section
        st.markdown("**Add Comment:**")
        comment_content = st.text_area(
            "Your comment",
            placeholder="Share your thoughts...",
            key=f"comment_content_{thread_id}",
            height=80
        )
        
        if st.button("💬 Add Comment", key=f"add_comment_{thread_id}"):
            if comment_content.strip():
                _add_comment(thread_id, comment_content.strip(), api_base)
            else:
                st.error("Please enter a comment.")


def _render_comment(comment: Dict, thread_id: str, api_base: str, level: int = 0) -> None:
    """Render an individual comment with optional nested replies."""
    # Indentation for nested comments
    indent = "　" * level  # Using full-width space for indentation
    
    # Comment content
    content = comment.get("content", "")
    created_at = comment.get("created_at", "")
    
    if created_at:
        try:
            dt = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
            time_str = dt.strftime("%m/%d %H:%M")
        except (ValueError, AttributeError):
            time_str = created_at[:16]
    else:
        time_str = "Unknown"
    
    # Display comment
    st.markdown(f"{indent}💬 **{time_str}:** {content}")
    
    # Reply button (only for top-level comments to keep UI simple)
    if level == 0:
        reply_key = f"reply_to_{comment['id']}"
        if st.button(f"{indent}↳ Reply", key=reply_key):
            st.session_state[f"show_reply_{comment['id']}"] = True
        
        # Reply form
        if st.session_state.get(f"show_reply_{comment['id']}", False):
            reply_content = st.text_area(
                "Reply",
                placeholder="Your reply...",
                key=f"reply_content_{comment['id']}",
                height=60
            )
            
            col1, col2 = st.columns([1, 1])
            with col1:
                if st.button("Send Reply", key=f"send_reply_{comment['id']}"):
                    if reply_content.strip():
                        _add_comment(thread_id, reply_content.strip(), api_base, parent_id=comment["id"])
                        st.session_state[f"show_reply_{comment['id']}"] = False
                        st.rerun()
                    else:
                        st.error("Please enter a reply.")
            
            with col2:
                if st.button("Cancel", key=f"cancel_reply_{comment['id']}"):
                    st.session_state[f"show_reply_{comment['id']}"] = False
                    st.rerun()
    
    # Render nested replies (if any)
    if comment.get("replies"):
        for reply in comment["replies"]:
            _render_comment(reply, thread_id, api_base, level + 1)


def _update_thread_status(thread_id: str, status: str, api_base: str) -> None:
    """Update thread status via API."""
    try:
        with st.spinner(f"Updating thread to {status}..."):
            response = httpx.patch(
                f"{api_base}/api/v1/threads/{thread_id}",
                json={"status": status}
            )
            response.raise_for_status()
        
        st.success(f"✅ Thread marked as {status}")
        st.rerun()
    except httpx.ConnectError:
        st.warning("⚠️ API server unavailable. Status update disabled.")
    except Exception as e:
        st.error(f"Failed to update thread: {e}")


def _add_comment(thread_id: str, content: str, api_base: str, parent_id: Optional[str] = None) -> None:
    """Add a comment to a thread via API."""
    try:
        with st.spinner("Adding comment..."):
            payload = {"content": content}
            if parent_id:
                payload["parent_id"] = parent_id
            
            response = httpx.post(
                f"{api_base}/api/v1/threads/{thread_id}/comments",
                json=payload
            )
            response.raise_for_status()
        
        st.success("✅ Comment added successfully!")
        st.rerun()
    except httpx.ConnectError:
        st.warning("⚠️ API server unavailable. Comment disabled.")
    except Exception as e:
        st.error(f"Failed to add comment: {e}")


def _get_demo_threads() -> List[Dict]:
    """Return demo threads data when API is unavailable."""
    return [
        {
            "id": str(uuid4()),
            "title": "Question about data preprocessing",
            "cell_index": 2,
            "status": "open",
            "created_at": "2025-01-01T10:00:00Z",
            "comments": [
                {
                    "id": str(uuid4()),
                    "content": "Could you explain the normalization method used here?",
                    "created_at": "2025-01-01T10:05:00Z",
                    "replies": []
                },
                {
                    "id": str(uuid4()),
                    "content": "I'm using z-score normalization to ensure all features have mean=0 and std=1.",
                    "created_at": "2025-01-01T10:15:00Z",
                    "replies": []
                }
            ]
        },
        {
            "id": str(uuid4()),
            "title": "Overall methodology review",
            "cell_index": None,
            "status": "resolved",
            "created_at": "2025-01-01T09:00:00Z",
            "comments": [
                {
                    "id": str(uuid4()),
                    "content": "The approach looks solid. Just minor formatting suggestions.",
                    "created_at": "2025-01-01T09:05:00Z",
                    "replies": []
                }
            ]
        }
    ]
