"""
Collapsible grouped sidebar navigation component.

Provides organized navigation with expandable groups.
"""

from __future__ import annotations

from typing import Optional

import streamlit as st

from scripts.dashboard.core.config import PAGE_GROUPS, GROUP_ORDER, GROUP_ICONS, PAGE_REGISTRY


def render_grouped_sidebar(current_page: Optional[str] = None) -> Optional[str]:
    """
    Render grouped sidebar navigation with collapsible sections.

    Args:
        current_page: Currently active page name for highlighting

    Returns:
        Selected page name if navigation occurred, None otherwise
    """
    selected_page = None
    
    for group in GROUP_ORDER:
        pages = PAGE_GROUPS.get(group, [])
        
        # Skip empty "Other" group
        if not pages and group == "Other":
            continue
        
        # Get group icon
        icon = GROUP_ICONS.get(group, "📁")
        
        # Get expansion state from session
        expansion_key = f"nav_group_{group}_expanded"
        is_expanded = st.session_state.get(expansion_key, group == "Home")
        
        # Render expandable group
        with st.expander(f"{icon} {group} ({len(pages)})", expanded=is_expanded):
            for page in pages:
                # Check if this page exists in registry
                if page not in PAGE_REGISTRY:
                    continue
                
                # Highlight current page
                button_type = "primary" if page == current_page else "secondary"
                
                # Navigation button
                button_key = f"nav_button_{group}_{page}"
                if st.button(
                    page,
                    key=button_key,
                    use_container_width=True,
                    type=button_type if page == current_page else None,
                ):
                    selected_page = page
                    # Update expansion state
                    st.session_state[expansion_key] = True
    
    return selected_page


def get_page_group(page_name: str) -> Optional[str]:
    """
    Get the group name for a given page.

    Args:
        page_name: Name of the page

    Returns:
        Group name if found, None otherwise
    """
    for group, pages in PAGE_GROUPS.items():
        if page_name in pages:
            return group
    return None


def get_group_page_count(group: str) -> int:
    """
    Get the number of pages in a group.

    Args:
        group: Group name

    Returns:
        Number of pages in the group
    """
    pages = PAGE_GROUPS.get(group, [])
    # Only count pages that exist in registry
    return len([p for p in pages if p in PAGE_REGISTRY])

