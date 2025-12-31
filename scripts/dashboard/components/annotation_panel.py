"""Annotation panel component for displaying and managing inline annotations."""

from __future__ import annotations

import requests
import streamlit as st
from datetime import datetime
from typing import Any, Dict, Optional
from uuid import UUID

from amprenta_rag.auth.session import get_current_user


def render_annotation_panel(
    entity_type: str, 
    entity_id: str | UUID,
    position_type: Optional[str] = None,
    position_data: Optional[Dict[str, Any]] = None,
    show_create_form: bool = True
) -> None:
    """
    Render an annotation panel for an entity or specific position.

    Args:
        entity_type: Type of entity ("notebook", "dataset", "experiment")
        entity_id: UUID string or UUID object of the entity
        position_type: Optional position type filter
        position_data: Optional position data for new annotations
        show_create_form: Whether to show the create annotation form
    """
    # Convert entity_id to UUID if it's a string
    if isinstance(entity_id, str):
        entity_id = UUID(entity_id)

    user = get_current_user()
    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None

    st.markdown("### 📝 Inline Annotations")

    # Filter controls
    col1, col2 = st.columns(2)
    with col1:
        status_filter = st.selectbox(
            "Status",
            options=["all", "open", "resolved"],
            index=0,
            key=f"annotation_status_filter_{entity_id}"
        )
    
    with col2:
        position_filter = st.selectbox(
            "Position Type",
            options=["all", "cell", "column", "row", "field", "range"],
            index=0,
            key=f"annotation_position_filter_{entity_id}"
        )

    # Get annotations from API
    annotations = _fetch_annotations(
        entity_type=entity_type,
        entity_id=entity_id,
        status_filter=None if status_filter == "all" else status_filter,
        position_type_filter=None if position_filter == "all" else position_filter
    )

    if annotations:
        annotation_list = annotations.get("annotations", [])
        counts = {
            "total": annotations.get("total", 0),
            "open": annotations.get("open_count", 0),
            "resolved": annotations.get("resolved_count", 0)
        }
        
        # Display counts
        st.metric("Total Annotations", counts["total"])
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Open", counts["open"])
        with col2:
            st.metric("Resolved", counts["resolved"])

        st.markdown("---")

        # Display annotations
        if annotation_list:
            for annotation in annotation_list:
                _render_annotation_card(annotation, user_id, entity_type, entity_id)
        else:
            st.info("No annotations match the current filters.")
    else:
        st.info("No annotations found for this entity.")

    # Create new annotation form
    if show_create_form and user_id:
        st.markdown("---")
        _render_create_annotation_form(
            entity_type=entity_type,
            entity_id=entity_id,
            position_type=position_type,
            position_data=position_data,
            user_id=user_id
        )
    elif show_create_form:
        st.info("Please log in to create annotations.")


def _fetch_annotations(
    entity_type: str,
    entity_id: UUID,
    status_filter: Optional[str] = None,
    position_type_filter: Optional[str] = None
) -> Optional[Dict]:
    """Fetch annotations from the API."""
    try:
        # Build query parameters
        params = {
            "entity_type": entity_type,
            "entity_id": str(entity_id)
        }
        
        if status_filter:
            params["status"] = status_filter
        if position_type_filter:
            params["position_type"] = position_type_filter

        # Make API request
        response = requests.get(
            "http://localhost:8000/api/v1/annotations",
            params=params,
            timeout=10
        )
        
        if response.status_code == 200:
            return response.json()
        else:
            st.error(f"Failed to fetch annotations: {response.status_code}")
            return None
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching annotations: {e}")
        return None


def _render_annotation_card(annotation: Dict, user_id: Optional[UUID], entity_type: str, entity_id: UUID) -> None:
    """Render a single annotation card."""
    annotation_id = annotation["id"]
    status = annotation["status"]
    
    # Status badge
    status_color = "🟢" if status == "open" else "🔵"
    status_text = f"{status_color} {status.title()}"
    
    with st.expander(f"{status_text} - {_format_position(annotation)}", expanded=(status == "open")):
        # Content
        st.markdown(f"**Content:** {annotation['content']}")
        
        # Metadata
        col1, col2 = st.columns(2)
        with col1:
            created_at = datetime.fromisoformat(annotation['created_at'].replace('Z', '+00:00'))
            st.caption(f"Created: {created_at.strftime('%Y-%m-%d %H:%M')}")
        
        with col2:
            if annotation.get('resolved_at'):
                resolved_at = datetime.fromisoformat(annotation['resolved_at'].replace('Z', '+00:00'))
                st.caption(f"Resolved: {resolved_at.strftime('%Y-%m-%d %H:%M')}")

        # Action buttons
        if user_id:
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if status == "open":
                    if st.button("✅ Resolve", key=f"resolve_{annotation_id}"):
                        _resolve_annotation(annotation_id)
                        st.rerun()
                else:
                    if st.button("🔄 Reopen", key=f"reopen_{annotation_id}"):
                        _reopen_annotation(annotation_id)
                        st.rerun()
            
            with col2:
                if st.button("💬 Reply", key=f"reply_{annotation_id}"):
                    st.session_state[f"show_reply_{annotation_id}"] = True
                    st.rerun()
            
            with col3:
                if st.button("🗑️ Delete", key=f"delete_{annotation_id}"):
                    _delete_annotation(annotation_id)
                    st.rerun()

        # Reply form
        if st.session_state.get(f"show_reply_{annotation_id}", False):
            _render_reply_form(annotation_id, user_id)

        # Show replies if any
        if annotation.get("parent_id") is None:
            _show_replies(annotation_id, entity_type, entity_id)


def _render_create_annotation_form(
    entity_type: str,
    entity_id: UUID,
    position_type: Optional[str] = None,
    position_data: Optional[Dict[str, Any]] = None,
    user_id: Optional[UUID] = None
) -> None:
    """Render the create annotation form."""
    st.markdown("#### ➕ Create New Annotation")
    
    with st.form(f"create_annotation_{entity_id}", clear_on_submit=True):
        # Position type selection
        pos_type = st.selectbox(
            "Position Type",
            options=["cell", "column", "row", "field", "range"],
            index=0 if position_type is None else ["cell", "column", "row", "field", "range"].index(position_type),
            disabled=position_type is not None
        )
        
        # Position data fields based on type
        pos_data = {}
        
        if pos_type == "cell":
            cell_index = st.number_input("Cell Index", min_value=0, value=position_data.get("cell_index", 0) if position_data else 0)
            pos_data = {"cell_index": int(cell_index)}
        elif pos_type == "column":
            column_name = st.text_input("Column Name", value=position_data.get("column", "") if position_data else "")
            pos_data = {"column": column_name}
        elif pos_type == "row":
            row_index = st.number_input("Row Index", min_value=0, value=position_data.get("row_index", 0) if position_data else 0)
            pos_data = {"row_index": int(row_index)}
        elif pos_type == "field":
            field_name = st.text_input("Field Name", value=position_data.get("field", "") if position_data else "")
            pos_data = {"field": field_name}
        elif pos_type == "range":
            col1, col2 = st.columns(2)
            with col1:
                start_cell = st.number_input("Start Cell", min_value=0, value=position_data.get("start_cell", 0) if position_data else 0)
            with col2:
                end_cell = st.number_input("End Cell", min_value=0, value=position_data.get("end_cell", 0) if position_data else 0)
            pos_data = {"start_cell": int(start_cell), "end_cell": int(end_cell)}
        
        # Content
        content = st.text_area("Annotation Content", placeholder="Describe the issue or observation...", height=100)
        
        # Submit button
        submitted = st.form_submit_button("📝 Create Annotation", type="primary")
        
        if submitted:
            if content.strip():
                success = _create_annotation(
                    entity_type=entity_type,
                    entity_id=entity_id,
                    position_type=pos_type,
                    position_data=pos_data,
                    content=content.strip()
                )
                
                if success:
                    st.success("Annotation created successfully!")
                    st.rerun()
            else:
                st.warning("Please enter annotation content.")


def _render_reply_form(annotation_id: str, user_id: Optional[UUID]) -> None:
    """Render reply form for an annotation."""
    with st.form(f"reply_form_{annotation_id}"):
        reply_content = st.text_area("Reply", placeholder="Write your reply...", height=80)
        
        col1, col2 = st.columns(2)
        with col1:
            submitted = st.form_submit_button("💬 Add Reply", type="primary")
        with col2:
            cancelled = st.form_submit_button("Cancel")
        
        if cancelled:
            st.session_state[f"show_reply_{annotation_id}"] = False
            st.rerun()
        
        if submitted:
            if reply_content.strip():
                success = _create_reply(annotation_id, reply_content.strip())
                if success:
                    st.success("Reply added!")
                    st.session_state[f"show_reply_{annotation_id}"] = False
                    st.rerun()
            else:
                st.warning("Please enter reply content.")


def _show_replies(annotation_id: str, entity_type: str, entity_id: UUID) -> None:
    """Show replies for an annotation."""
    # This would fetch replies from the API
    # For now, we'll show a placeholder
    pass


def _format_position(annotation: Dict) -> str:
    """Format position data for display."""
    position_type = annotation["position_type"]
    position_data = annotation["position_data"]
    
    if position_type == "cell":
        return f"Cell {position_data.get('cell_index', 'N/A')}"
    elif position_type == "column":
        return f"Column '{position_data.get('column', 'N/A')}'"
    elif position_type == "row":
        return f"Row {position_data.get('row_index', 'N/A')}"
    elif position_type == "field":
        return f"Field '{position_data.get('field', 'N/A')}'"
    elif position_type == "range":
        start = position_data.get('start_cell', 'N/A')
        end = position_data.get('end_cell', 'N/A')
        return f"Cells {start}-{end}"
    else:
        return f"{position_type}: {position_data}"


def _create_annotation(
    entity_type: str,
    entity_id: UUID,
    position_type: str,
    position_data: Dict[str, Any],
    content: str
) -> bool:
    """Create a new annotation via API."""
    try:
        payload = {
            "entity_type": entity_type,
            "entity_id": str(entity_id),
            "position_type": position_type,
            "position_data": position_data,
            "content": content,
            "parent_id": None
        }
        
        response = requests.post(
            "http://localhost:8000/api/v1/annotations",
            json=payload,
            timeout=10
        )
        
        if response.status_code == 201:
            return True
        else:
            st.error(f"Failed to create annotation: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error creating annotation: {e}")
        return False


def _resolve_annotation(annotation_id: str) -> bool:
    """Resolve an annotation via API."""
    try:
        response = requests.patch(
            f"http://localhost:8000/api/v1/annotations/{annotation_id}/resolve",
            timeout=10
        )
        
        if response.status_code == 200:
            st.success("Annotation resolved!")
            return True
        else:
            st.error(f"Failed to resolve annotation: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error resolving annotation: {e}")
        return False


def _reopen_annotation(annotation_id: str) -> bool:
    """Reopen an annotation via API."""
    try:
        response = requests.patch(
            f"http://localhost:8000/api/v1/annotations/{annotation_id}/reopen",
            timeout=10
        )
        
        if response.status_code == 200:
            st.success("Annotation reopened!")
            return True
        else:
            st.error(f"Failed to reopen annotation: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error reopening annotation: {e}")
        return False


def _delete_annotation(annotation_id: str) -> bool:
    """Delete an annotation via API."""
    try:
        response = requests.delete(
            f"http://localhost:8000/api/v1/annotations/{annotation_id}",
            timeout=10
        )
        
        if response.status_code == 204:
            st.success("Annotation deleted!")
            return True
        else:
            st.error(f"Failed to delete annotation: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error deleting annotation: {e}")
        return False


def _create_reply(annotation_id: str, content: str) -> bool:
    """Create a reply to an annotation via API."""
    try:
        payload = {"content": content}
        
        response = requests.post(
            f"http://localhost:8000/api/v1/annotations/{annotation_id}/reply",
            json=payload,
            timeout=10
        )
        
        if response.status_code == 201:
            return True
        else:
            st.error(f"Failed to create reply: {response.status_code}")
            return False
            
    except requests.exceptions.RequestException as e:
        st.error(f"Error creating reply: {e}")
        return False


def render_annotation_indicator(
    entity_type: str,
    entity_id: UUID,
    position_type: str,
    position_data: Dict[str, Any],
    label: str = "📝"
) -> bool:
    """
    Render an annotation indicator button that opens the annotation panel.
    
    Args:
        entity_type: Type of entity
        entity_id: Entity UUID
        position_type: Position type
        position_data: Position data
        label: Button label
        
    Returns:
        True if button was clicked
    """
    # Get annotation count for this position
    annotations = _fetch_annotations(
        entity_type=entity_type,
        entity_id=entity_id,
        position_type_filter=position_type
    )
    
    count = 0
    if annotations:
        # Filter by exact position match
        for ann in annotations.get("annotations", []):
            if ann["position_data"] == position_data:
                count += 1
    
    # Show button with count badge
    button_text = f"{label} ({count})" if count > 0 else label
    button_key = f"annotation_indicator_{entity_type}_{entity_id}_{position_type}_{hash(str(position_data))}"
    
    if st.button(button_text, key=button_key, help="View/Add annotations"):
        # Store the position context for the annotation panel
        st.session_state["annotation_context"] = {
            "entity_type": entity_type,
            "entity_id": entity_id,
            "position_type": position_type,
            "position_data": position_data
        }
        return True
    
    return False
