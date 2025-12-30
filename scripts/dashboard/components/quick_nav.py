"""Quick navigation command palette component."""

from __future__ import annotations

import streamlit as st
from scripts.dashboard.core.config import PAGE_REGISTRY, PAGE_GROUPS, GROUP_ICONS


def fuzzy_match(query: str, page_name: str) -> float:
    """
    Simple fuzzy match score (0-1). Higher is better.
    
    Args:
        query: Search query string
        page_name: Page name to match against
        
    Returns:
        Match score between 0 and 1
    """
    if not query or not page_name:
        return 0.0
    
    query = query.lower().strip()
    page = page_name.lower()
    
    # Exact substring match gets high score
    if query in page:
        # Perfect match
        if query == page:
            return 1.0
        # Substring match - score based on how much of the page name matches
        return 0.8 + (len(query) / len(page)) * 0.2
    
    # Check word starts for acronym-style matching
    words = page.split()
    if len(words) > 1:
        # Check if query matches first letters of words (e.g., "gc" -> "Generative Chemistry")
        first_letters = ''.join(w[0] for w in words if w).lower()
        if query in first_letters:
            return 0.6
        
        # Check if any word starts with the query
        word_matches = sum(1 for w in words if w.lower().startswith(query))
        if word_matches > 0:
            return 0.4 + (word_matches / len(words)) * 0.2
    
    # Check character-by-character fuzzy matching
    query_chars = list(query)
    page_chars = list(page)
    matches = 0
    page_idx = 0
    
    for q_char in query_chars:
        while page_idx < len(page_chars):
            if page_chars[page_idx] == q_char:
                matches += 1
                page_idx += 1
                break
            page_idx += 1
        else:
            break  # Couldn't find this character
    
    if matches == len(query_chars):
        return 0.3 * (matches / len(page_chars))
    
    return 0.0


def search_pages(query: str, limit: int = 10) -> list[dict]:
    """
    Search pages by fuzzy match.
    
    Args:
        query: Search query string
        limit: Maximum number of results to return
        
    Returns:
        List of page dictionaries with metadata
    """
    if not query:
        return []
    
    results = []
    for page in PAGE_REGISTRY.keys():
        score = fuzzy_match(query, page)
        if score > 0.1:  # Only include decent matches
            group = None
            for g, pages in PAGE_GROUPS.items():
                if page in pages:
                    group = g
                    break
            
            results.append({
                "page": page,
                "group": group,
                "icon": GROUP_ICONS.get(group, "📄"),
                "score": score,
            })
    
    # Sort by score (highest first)
    results.sort(key=lambda x: x["score"], reverse=True)
    return results[:limit]


def render_command_palette() -> None:
    """Render command palette modal (triggered by session state)."""
    if not st.session_state.get("show_command_palette"):
        return
    
    # Create a modal-like container
    with st.container():
        st.markdown("---")
        st.markdown("### 🔍 Quick Navigation")
        st.caption("Type to search pages, press Enter to select first result")
        
        # Search input
        query = st.text_input(
            "Search pages...",
            key="command_palette_query",
            placeholder="e.g., 'chemistry', 'gc' (Generative Chemistry), 'admin'",
            help="Use Ctrl+K to open, Esc to close",
        )
        
        # Auto-focus on search input when opened
        if query is None:
            st.session_state["command_palette_query"] = ""
        
        # Search results
        results = search_pages(query)
        
        if results:
            st.markdown("**Results:**")
            
            # Limit display to prevent overwhelming UI
            display_results = results[:8]
            
            for i, r in enumerate(display_results):
                # Create display text with icon and group
                display = f"{r['icon']} {r['page']}"
                if r['group']:
                    display += f" • {r['group']}"
                
                # Add score for debugging (can be removed in production)
                score_text = f" ({r['score']:.2f})" if st.session_state.get("debug_scores") else ""
                
                # First result gets primary styling for Enter key selection
                button_type = "primary" if i == 0 else "secondary"
                
                if st.button(
                    display + score_text, 
                    key=f"qnav_{r['page']}", 
                    use_container_width=True,
                    type=button_type
                ):
                    st.session_state["selected_page"] = r["page"]
                    st.session_state["show_command_palette"] = False
                    st.session_state["command_palette_query"] = ""  # Clear search
                    st.rerun()
            
            if len(results) > len(display_results):
                st.caption(f"... and {len(results) - len(display_results)} more results")
                
        elif query:
            st.info(f"No pages found for '{query}'")
            st.caption("Try shorter keywords or check spelling")
        else:
            # Show recent pages or suggestions when no query
            recent_pages = st.session_state.get("recent_pages", [])
            if recent_pages:
                st.markdown("**Recent Pages:**")
                for page in recent_pages[:5]:
                    if st.button(f"🕐 {page}", key=f"recent_qnav_{page}", use_container_width=True):
                        st.session_state["selected_page"] = page
                        st.session_state["show_command_palette"] = False
                        st.session_state["command_palette_query"] = ""
                        st.rerun()
        
        # Control buttons
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col1:
            if st.button("Close", key="close_command_palette"):
                st.session_state["show_command_palette"] = False
                st.session_state["command_palette_query"] = ""
                st.rerun()
        
        with col2:
            # Debug toggle for development (only in debug mode)
            if os.environ.get("AMPRENTA_DEBUG"):
                if st.button("Debug", key="toggle_debug"):
                    st.session_state["debug_scores"] = not st.session_state.get("debug_scores", False)
                    st.rerun()
        
        with col3:
            # Clear search
            if st.button("Clear", key="clear_search"):
                st.session_state["command_palette_query"] = ""
                st.rerun()
        
        st.markdown("---")


def render_pinned_pages() -> None:
    """Render pinned quick-access pages at top of sidebar."""
    pinned = st.session_state.get("pinned_pages", [])
    
    if not pinned:
        # Show default pinned pages for first-time users
        default_pinned = ["Overview", "Chemistry", "Experiments", "Analysis Tools"]
        pinned = [p for p in default_pinned if p in PAGE_REGISTRY]
        if pinned:
            st.session_state["pinned_pages"] = pinned
    
    if not pinned:
        return
    
    st.markdown("#### 📌 Quick Access")
    
    # Create columns for horizontal layout if there are few pinned pages
    if len(pinned) <= 2:
        cols = st.columns(len(pinned))
        for i, page in enumerate(pinned):
            if page in PAGE_REGISTRY:
                with cols[i]:
                    if st.button(page, key=f"pinned_{page}", use_container_width=True):
                        st.session_state["selected_page"] = page
                        st.rerun()
    else:
        # Vertical layout for many pinned pages
        for page in pinned:
            if page in PAGE_REGISTRY:
                if st.button(page, key=f"pinned_{page}", use_container_width=True):
                    st.session_state["selected_page"] = page
                    st.rerun()
    
    # Pin management
    with st.expander("⚙️ Manage Pins", expanded=False):
        st.caption("Pin frequently used pages for quick access")
        
        # Add new pin
        available_pages = [p for p in PAGE_REGISTRY.keys() if p not in pinned]
        if available_pages:
            new_pin = st.selectbox(
                "Add page to pins:",
                [""] + sorted(available_pages),
                key="new_pin_select"
            )
            if st.button("📌 Pin", key="add_pin") and new_pin:
                pinned.append(new_pin)
                st.session_state["pinned_pages"] = pinned
                st.rerun()
        
        # Remove pins
        if pinned:
            remove_pin = st.selectbox(
                "Remove pin:",
                [""] + pinned,
                key="remove_pin_select"
            )
            if st.button("🗑️ Remove", key="remove_pin") and remove_pin:
                pinned.remove(remove_pin)
                st.session_state["pinned_pages"] = pinned
                st.rerun()
    
    st.divider()


def inject_quick_nav_js() -> None:
    """Inject JavaScript for Ctrl+K keyboard shortcut."""
    js_code = """
    <script>
    document.addEventListener('keydown', function(e) {
        // Ctrl+K or Cmd+K to open command palette
        if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
            e.preventDefault();
            
            // Set session state to show command palette
            // This is a workaround since we can't directly modify Streamlit session state from JS
            // In practice, this would need to be handled differently
            console.log('Quick nav shortcut triggered');
            
            // Try to focus search input if it exists
            const searchInput = document.querySelector('input[placeholder*="Search pages"]');
            if (searchInput) {
                searchInput.focus();
            }
        }
        
        // Escape to close command palette
        if (e.key === 'Escape') {
            // This would close the command palette
            console.log('Escape key pressed');
        }
    });
    </script>
    """
    st.markdown(js_code, unsafe_allow_html=True)


def toggle_command_palette() -> None:
    """Toggle command palette visibility."""
    current_state = st.session_state.get("show_command_palette", False)
    st.session_state["show_command_palette"] = not current_state
    
    if st.session_state["show_command_palette"]:
        # Clear previous search when opening
        st.session_state["command_palette_query"] = ""


def get_page_suggestions(current_page: str) -> list[str]:
    """
    Get page suggestions based on current page context.
    
    Args:
        current_page: Currently active page
        
    Returns:
        List of suggested page names
    """
    suggestions = []
    
    # Find current page's group
    current_group = None
    for group, pages in PAGE_GROUPS.items():
        if current_page in pages:
            current_group = group
            break
    
    # Suggest related pages from same group
    if current_group:
        group_pages = PAGE_GROUPS[current_group]
        suggestions.extend([p for p in group_pages if p != current_page][:3])
    
    # Add some popular pages
    popular_pages = ["Overview", "Chemistry", "Experiments", "Analysis Tools", "Company Settings"]
    for page in popular_pages:
        if page not in suggestions and page != current_page and page in PAGE_REGISTRY:
            suggestions.append(page)
    
    return suggestions[:5]
