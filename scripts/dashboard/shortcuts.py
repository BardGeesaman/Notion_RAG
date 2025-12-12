"""Keyboard shortcuts for the Streamlit dashboard."""
from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

SHORTCUTS = {
    "?": "Show keyboard shortcuts help",
    "g+h": "Go to Home/Overview",
    "g+e": "Go to Experiments",
    "g+c": "Go to Chemistry",
    "g+a": "Go to Analysis Tools",
    "g+s": "Go to Search",
    "g+d": "Go to Data Ingestion",
    "g+p": "Go to Protocols",
    "g+q": "Go to Q&A Tracker",
    "/": "Focus search",
    "esc": "Close modals/dialogs",
}


def render_shortcuts_help() -> None:
    """Display keyboard shortcuts help in an expander."""
    with st.expander("⌨️ Keyboard Shortcuts", expanded=False):
        st.markdown("**Available shortcuts:**")
        for shortcut, description in SHORTCUTS.items():
            st.markdown(f"- `{shortcut}` - {description}")


def inject_shortcuts_js() -> None:
    """
    Inject JavaScript to handle keyboard shortcuts.
    
    Uses st.components.v1.html to inject JavaScript that listens for keydown events
    and updates Streamlit session state for navigation.
    """
    js_code = """
    <script>
    (function() {
        let keys = [];
        let lastKeyTime = Date.now();
        
        function handleKeyDown(event) {
            const now = Date.now();
            
            // Reset keys if too much time has passed (for multi-key shortcuts)
            if (now - lastKeyTime > 1000) {
                keys = [];
            }
            lastKeyTime = now;
            
            // Don't trigger shortcuts when typing in inputs
            if (event.target.tagName === 'INPUT' || 
                event.target.tagName === 'TEXTAREA' || 
                event.target.isContentEditable) {
                return;
            }
            
            const key = event.key.toLowerCase();
            
            // Single key shortcuts
            if (key === '?' && !event.shiftKey && !event.ctrlKey && !event.metaKey) {
                event.preventDefault();
                window.parent.postMessage({
                    type: 'streamlit:setComponentValue',
                    value: {action: 'show_shortcuts'}
                }, '*');
                return;
            }
            
            if (key === '/' && !event.shiftKey && !event.ctrlKey && !event.metaKey) {
                event.preventDefault();
                window.parent.postMessage({
                    type: 'streamlit:setComponentValue',
                    value: {action: 'focus_search'}
                }, '*');
                return;
            }
            
            if (key === 'escape' || key === 'esc') {
                event.preventDefault();
                window.parent.postMessage({
                    type: 'streamlit:setComponentValue',
                    value: {action: 'close_modals'}
                }, '*');
                return;
            }
            
            // Multi-key shortcuts (g+key pattern)
            if (key === 'g') {
                keys.push('g');
                return;
            }
            
            if (keys.length > 0 && keys[0] === 'g') {
                let action = null;
                
                if (key === 'h') {
                    action = 'navigate_home';
                } else if (key === 'e') {
                    action = 'navigate_experiments';
                } else if (key === 'c') {
                    action = 'navigate_chemistry';
                } else if (key === 'a') {
                    action = 'navigate_analysis';
                } else if (key === 's') {
                    action = 'navigate_search';
                } else if (key === 'd') {
                    action = 'navigate_ingestion';
                } else if (key === 'p') {
                    action = 'navigate_protocols';
                } else if (key === 'q') {
                    action = 'navigate_qa';
                }
                
                if (action) {
                    event.preventDefault();
                    window.parent.postMessage({
                        type: 'streamlit:setComponentValue',
                        value: {action: action}
                    }, '*');
                }
                
                keys = [];
            } else {
                keys = [];
            }
        }
        
        // Listen for keydown events
        document.addEventListener('keydown', handleKeyDown);
        
        // Cleanup function (not really needed for Streamlit but good practice)
        window.addEventListener('beforeunload', function() {
            document.removeEventListener('keydown', handleKeyDown);
        });
    })();
    </script>
    """
    
    components.html(js_code, height=0)


def handle_shortcut_action(action: str) -> None:
    """
    Handle shortcut actions by updating session state.
    
    Args:
        action: Action string from JavaScript
    """
    page_mapping = {
        'navigate_home': 'Overview',
        'navigate_experiments': 'Experiments',
        'navigate_chemistry': 'Chemistry',
        'navigate_analysis': 'Analysis Tools',
        'navigate_search': 'Search',
        'navigate_ingestion': 'Data Ingestion',
        'navigate_protocols': 'Protocols',
        'navigate_qa': 'Q&A Tracker',
    }
    
    if action in page_mapping:
        st.session_state["selected_page"] = page_mapping[action]
        st.rerun()
    elif action == 'show_shortcuts':
        st.session_state["show_shortcuts"] = True
    elif action == 'focus_search':
        st.session_state["focus_search"] = True
    elif action == 'close_modals':
        # Close any open modals/expanders
        if "show_shortcuts" in st.session_state:
            st.session_state.pop("show_shortcuts", None)
        if "show_register" in st.session_state:
            st.session_state.pop("show_register", None)
        st.rerun()
