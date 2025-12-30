"""Breadcrumb navigation component."""

from __future__ import annotations

import streamlit as st
from scripts.dashboard.core.config import PAGE_GROUPS, GROUP_ICONS
from scripts.dashboard.components.sidebar_nav import get_page_group


def get_breadcrumb_path(current_page: str) -> list[tuple[str, str, str | None]]:
    """
    Get breadcrumb path for current page.
    
    Args:
        current_page: Name of the current page
    
    Returns:
        List of (icon, label, page_name) tuples
        e.g., [("🏠", "Home", "Overview"), ("⚗️", "Chemistry", None), ("", "Generative Chemistry", "Generative Chemistry")]
        
        Where:
        - icon: Emoji icon for the breadcrumb item
        - label: Display text for the breadcrumb item
        - page_name: Page to navigate to (None if not clickable)
    """
    # Special case for Home group pages - show minimal breadcrumb
    home_pages = PAGE_GROUPS.get("Home", [])
    if current_page in home_pages:
        if current_page == "Overview":
            return [("🏠", "Home", None), ("", "Overview", "Overview")]
        else:
            return [("🏠", "Home", "Overview"), ("", current_page, current_page)]
    
    # Find which group contains this page
    for group, pages in PAGE_GROUPS.items():
        if current_page in pages:
            icon = GROUP_ICONS.get(group, "📁")
            return [
                ("🏠", "Home", "Overview"),
                (icon, group, None),  # Group not clickable
                ("", current_page, current_page),
            ]
    
    # Fallback if page not in any group
    return [("🏠", "Home", "Overview"), ("", current_page, current_page)]


def render_breadcrumbs(current_page: str) -> None:
    """
    Render breadcrumb navigation bar.
    
    Args:
        current_page: Name of the current page to show breadcrumbs for
    """
    if not current_page:
        return
    
    path = get_breadcrumb_path(current_page)
    
    if len(path) <= 1:
        return  # Don't show breadcrumbs for single-level paths
    
    # Calculate number of columns needed (items + separators)
    num_cols = len(path) * 2 - 1
    cols = st.columns(num_cols)
    col_idx = 0
    
    for i, (icon, label, page_name) in enumerate(path):
        with cols[col_idx]:
            display = f"{icon} {label}" if icon else label
            
            if page_name and page_name != current_page:
                # Clickable link for navigation
                if st.button(
                    display, 
                    key=f"breadcrumb_{i}_{page_name}", 
                    type="secondary",
                    use_container_width=True
                ):
                    st.session_state["selected_page"] = page_name
                    st.rerun()
            else:
                # Current page (not clickable) or group label
                if page_name == current_page:
                    # Current page - style differently
                    st.markdown(f"**{display}**")
                else:
                    # Group label - muted style
                    st.markdown(f"*{display}*")
        
        col_idx += 1
        
        # Add separator except after last item
        if i < len(path) - 1 and col_idx < len(cols):
            with cols[col_idx]:
                st.markdown("›", help="Breadcrumb separator")
            col_idx += 1




def get_group_icon(group_name: str) -> str:
    """
    Get the icon for a given group.
    
    Args:
        group_name: Name of the group
        
    Returns:
        Icon emoji for the group
    """
    return GROUP_ICONS.get(group_name, "📁")


def render_compact_breadcrumbs(current_page: str) -> None:
    """
    Render compact breadcrumb navigation in a single line.
    
    Args:
        current_page: Name of the current page to show breadcrumbs for
    """
    if not current_page:
        return
    
    path = get_breadcrumb_path(current_page)
    
    if len(path) <= 1:
        return
    
    # Build breadcrumb string
    breadcrumb_parts = []
    
    for i, (icon, label, page_name) in enumerate(path):
        display = f"{icon} {label}" if icon else label
        
        if page_name and page_name != current_page:
            # Make it a clickable link (using markdown link syntax won't work in Streamlit)
            # So we'll just show it as text for compact mode
            breadcrumb_parts.append(display)
        else:
            if page_name == current_page:
                # Current page - bold
                breadcrumb_parts.append(f"**{display}**")
            else:
                # Group label - italic
                breadcrumb_parts.append(f"*{display}*")
    
    # Join with separators
    breadcrumb_text = " › ".join(breadcrumb_parts)
    st.markdown(breadcrumb_text)
