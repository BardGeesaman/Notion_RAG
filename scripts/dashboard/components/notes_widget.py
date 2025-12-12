"""Notes widget component for embedding in detail pages."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.utils.notes import add_note, get_notes, delete_note
from scripts.dashboard.db_session import db_session
from amprenta_rag.auth.session import get_current_user


def render_notes_widget(entity_type: str, entity_id: str) -> None:
    """
    Render a notes widget for an entity.
    
    Args:
        entity_type: Type of entity (experiment, compound, signature, etc.)
        entity_id: UUID of the entity
    """
    user = get_current_user()
    
    st.markdown("### 📝 Notes")
    
    with db_session() as db:
        # Display existing notes
        notes = get_notes(entity_type, entity_id, db)
        
        if notes:
            for note in notes:
                with st.container():
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        author = note.created_by.username if note.created_by else "Unknown"
                        timestamp = note.created_at.strftime("%Y-%m-%d %H:%M") if note.created_at else "Unknown"
                        st.markdown(f"**{author}** ({timestamp})")
                        st.markdown(note.content)
                    with col2:
                        # Only allow deletion by the author or admin
                        if user and (user.get("id") == str(note.created_by_id) or user.get("role") == "admin"):
                            if st.button("🗑️", key=f"delete_note_{note.id}", help="Delete note"):
                                delete_note(str(note.id), db)
                                st.rerun()
                    st.divider()
        else:
            st.info("No notes yet. Add one below.")
        
        # Add new note
        if user:
            with st.form(f"add_note_{entity_id}", clear_on_submit=True):
                new_note = st.text_area("Add a note", placeholder="Enter your note here...", key=f"note_input_{entity_id}")
                if st.form_submit_button("Add Note", type="primary"):
                    if new_note.strip():
                        try:
                            add_note(entity_type, entity_id, new_note.strip(), user.get("id"), db)
                            st.success("Note added!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to add note: {e}")
                    else:
                        st.warning("Note cannot be empty")
        else:
            st.info("Please log in to add notes.")
