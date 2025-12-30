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
    current_group = get_page_group(current_page) if current_page else None
    
    for group in GROUP_ORDER:
        pages = PAGE_GROUPS.get(group, [])
        
        # Skip empty "Other" group
        if not pages and group == "Other":
            continue
        
        # Count valid pages (only those in registry)
        valid_pages = [p for p in pages if p in PAGE_REGISTRY]
        page_count = len(valid_pages)
        
        if page_count == 0:
            continue
        
        # Get group icon
        icon = GROUP_ICONS.get(group, "📁")
        
        # Get expansion state from session - keep current group expanded
        expansion_key = f"nav_group_{group}_expanded"
        default_expanded = group == "Home" or group == current_group
        is_expanded = st.session_state.get(expansion_key, default_expanded)
        
        # Check if group is pinned
        pin_key = f"nav_group_{group}_pinned"
        is_pinned = st.session_state.get(pin_key, False)
        
        # Pin icon for pinned groups
        pin_icon = "⭐" if is_pinned else ""
        
        # Highlight current page's group with different styling
        group_label = f"{icon} {group} ({page_count}){pin_icon}"
        
        # Store expansion state when expander is used
        with st.expander(group_label, expanded=is_expanded or is_pinned) as expander:
            # Update expansion state
            if expander:
                st.session_state[expansion_key] = True
            
            # Add pin/unpin button
            col1, col2 = st.columns([3, 1])
            with col2:
                pin_label = "📌" if not is_pinned else "📍"
                if st.button(pin_label, key=f"pin_{group}", help="Pin group"):
                    st.session_state[pin_key] = not is_pinned
                    st.rerun()
            
            # Render pages in this group
            for page in valid_pages:
                # Highlight current page
                button_type = "primary" if page == current_page else "secondary"
                
                # Add indicator for current page
                page_label = f"→ {page}" if page == current_page else page
                
                # Navigation button
                button_key = f"nav_button_{group}_{page}"
                if st.button(
                    page_label,
                    key=button_key,
                    use_container_width=True,
                    type=button_type,
                ):
                    selected_page = page
                    # Keep this group expanded after navigation
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

