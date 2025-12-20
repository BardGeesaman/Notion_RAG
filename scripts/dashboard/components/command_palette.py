"""Command palette component for quick actions."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.utils.actions import QUICK_ACTIONS, search_actions


def render_command_palette() -> None:
    """Render the command palette modal."""
    # Initialize palette state
    if "show_palette" not in st.session_state:
        st.session_state["show_palette"] = False

    # Show palette if requested
    if st.session_state.get("show_palette", False):
        # Create modal-like overlay
        with st.container():
            st.markdown("---")
            st.markdown("### ⌘ Command Palette")
            st.caption("Press Ctrl+K to open, Esc to close")

            # Search input
            search_query = st.text_input(
                "Search actions...",
                value="",
                key="palette_search",
                placeholder="Type to search...",
                autofocus=True,
            )  # type: ignore[call-overload]

            # Filter actions
            if search_query:
                actions = search_actions(search_query)
            else:
                actions = QUICK_ACTIONS

            # Display actions
            if actions:
                st.markdown(f"**{len(actions)} action(s) found**")
                for action in actions[:10]:  # Limit to 10 for display
                    if st.button(
                        f"{action['icon']} **{action['name']}** - {action['description']}",
                        key=f"palette_action_{action['name']}",
                        use_container_width=True,
                    ):
                        st.session_state["selected_page"] = action["page"]
                        st.session_state["show_palette"] = False
                        st.rerun()
            else:
                st.info("No actions found matching your search.")

            # Close button
            if st.button("Close (Esc)", key="close_palette", use_container_width=True):
                st.session_state["show_palette"] = False
                st.rerun()
