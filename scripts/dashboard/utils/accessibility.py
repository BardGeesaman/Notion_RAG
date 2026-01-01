"""Accessibility utilities for dashboard components.

This module provides utilities to improve WCAG AA compliance within Streamlit's limitations.
Focus is on critical user paths: authentication, navigation, and core functionality.
"""

import streamlit as st
from typing import Optional, Any
import uuid


def add_aria_label(element_key: str, label: str):
    """Add ARIA label to Streamlit element via custom HTML/JS injection.
    
    Args:
        element_key: Streamlit element key or test ID
        label: ARIA label text for screen readers
    """
    # Use a unique ID to avoid script conflicts
    script_id = f"aria_{uuid.uuid4().hex[:8]}"
    
    st.markdown(
        f"""
        <script id="{script_id}">
            (function() {{
                // Try multiple selectors to find the element
                const selectors = [
                    '[data-testid="{element_key}"]',
                    '[key="{element_key}"]',
                    '#{element_key}'
                ];
                
                let element = null;
                for (const selector of selectors) {{
                    element = document.querySelector(selector);
                    if (element) break;
                }}
                
                if (element) {{
                    element.setAttribute('aria-label', '{label}');
                    element.setAttribute('title', '{label}');
                }} else {{
                    // Retry after a short delay for dynamically rendered elements
                    setTimeout(() => {{
                        for (const selector of selectors) {{
                            element = document.querySelector(selector);
                            if (element) {{
                                element.setAttribute('aria-label', '{label}');
                                element.setAttribute('title', '{label}');
                                break;
                            }}
                        }}
                    }}, 100);
                }}
                
                // Clean up script
                document.getElementById('{script_id}').remove();
            }})();
        </script>
        """,
        unsafe_allow_html=True
    )


def accessible_button(
    label: str, 
    key: str, 
    aria_label: Optional[str] = None,
    help: Optional[str] = None,
    **kwargs
) -> bool:
    """Render button with explicit ARIA label and accessibility features.
    
    Args:
        label: Button text
        key: Unique key for the button
        aria_label: Custom ARIA label (defaults to label)
        help: Help text for additional context
        **kwargs: Additional arguments passed to st.button
    
    Returns:
        Boolean indicating if button was clicked
    """
    result = st.button(label, key=key, help=help, **kwargs)
    
    # Add ARIA label if provided or use label as fallback
    final_aria_label = aria_label or label
    add_aria_label(key, final_aria_label)
    
    return result


def accessible_selectbox(
    label: str,
    options: list,
    key: str,
    aria_label: Optional[str] = None,
    help: Optional[str] = None,
    **kwargs
) -> Any:
    """Render selectbox with ARIA support and accessibility features.
    
    Args:
        label: Selectbox label
        options: List of options
        key: Unique key for the selectbox
        aria_label: Custom ARIA label (defaults to label)
        help: Help text for additional context
        **kwargs: Additional arguments passed to st.selectbox
    
    Returns:
        Selected value
    """
    result = st.selectbox(label, options, key=key, help=help, **kwargs)
    
    # Add ARIA label if provided or use label as fallback
    final_aria_label = aria_label or f"Select {label.lower()}"
    add_aria_label(key, final_aria_label)
    
    return result


def accessible_text_input(
    label: str,
    key: str,
    aria_label: Optional[str] = None,
    aria_describedby: Optional[str] = None,
    help: Optional[str] = None,
    **kwargs
) -> str:
    """Render text input with accessibility features.
    
    Args:
        label: Input label
        key: Unique key for the input
        aria_label: Custom ARIA label (defaults to label)
        aria_describedby: ID of element describing the input
        help: Help text for additional context
        **kwargs: Additional arguments passed to st.text_input
    
    Returns:
        Input value
    """
    result = st.text_input(label, key=key, help=help, **kwargs)
    
    # Add ARIA label and description
    final_aria_label = aria_label or label
    add_aria_label(key, final_aria_label)
    
    if aria_describedby:
        # Add aria-describedby attribute
        script_id = f"describedby_{uuid.uuid4().hex[:8]}"
        st.markdown(
            f"""
            <script id="{script_id}">
                (function() {{
                    const element = document.querySelector('[data-testid="{key}"]') || 
                                  document.querySelector('[key="{key}"]');
                    if (element) {{
                        const input = element.querySelector('input');
                        if (input) {{
                            input.setAttribute('aria-describedby', '{aria_describedby}');
                        }}
                    }}
                    document.getElementById('{script_id}').remove();
                }})();
            </script>
            """,
            unsafe_allow_html=True
        )
    
    return result


def render_skip_link(target_id: str = "main-content"):
    """Render skip-to-content link for keyboard users.
    
    Args:
        target_id: ID of the main content area to skip to
    """
    st.markdown(
        f"""
        <style>
            .skip-link {{
                position: absolute;
                left: -9999px;
                z-index: 9999;
                background: #000;
                color: #fff;
                padding: 8px 16px;
                text-decoration: none;
                border-radius: 4px;
                font-weight: bold;
                border: 2px solid #fff;
                transition: left 0.2s ease;
            }}
            .skip-link:focus {{
                left: 10px !important;
                top: 10px;
            }}
            .skip-link:hover {{
                background: #333;
            }}
        </style>
        <a href="#{target_id}" class="skip-link" tabindex="1">
            Skip to main content
        </a>
        <div id="{target_id}" tabindex="-1"></div>
        """,
        unsafe_allow_html=True
    )


def announce_to_screen_reader(message: str, priority: str = "polite"):
    """Announce message to screen readers via live region.
    
    Args:
        message: Message to announce
        priority: 'polite' or 'assertive' for announcement priority
    """
    region_id = f"sr_announce_{uuid.uuid4().hex[:8]}"
    
    st.markdown(
        f"""
        <div id="{region_id}" 
             role="status" 
             aria-live="{priority}" 
             aria-atomic="true" 
             style="position: absolute; left: -9999px; width: 1px; height: 1px; overflow: hidden;">
            {message}
        </div>
        <script>
            // Remove announcement after 5 seconds to prevent accumulation
            setTimeout(() => {{
                const element = document.getElementById('{region_id}');
                if (element) element.remove();
            }}, 5000);
        </script>
        """,
        unsafe_allow_html=True
    )


def add_navigation_landmark(role: str = "navigation", label: Optional[str] = None):
    """Add navigation landmark for screen readers.
    
    Args:
        role: ARIA role (navigation, main, banner, etc.)
        label: ARIA label for the landmark
    """
    landmark_id = f"landmark_{uuid.uuid4().hex[:8]}"
    aria_label_attr = f'aria-label="{label}"' if label else ''
    
    st.markdown(
        f"""
        <div id="{landmark_id}" role="{role}" {aria_label_attr}>
        <!-- Navigation landmark start -->
        </div>
        """,
        unsafe_allow_html=True
    )


def add_heading_structure(text: str, level: int = 1, id: Optional[str] = None):
    """Add properly structured heading with accessibility features.
    
    Args:
        text: Heading text
        level: Heading level (1-6)
        id: Optional ID for the heading
    """
    if not 1 <= level <= 6:
        level = 1
    
    heading_id = id or f"heading_{uuid.uuid4().hex[:8]}"
    
    st.markdown(
        f"""
        <h{level} id="{heading_id}" tabindex="-1" style="margin-top: 1rem; margin-bottom: 0.5rem;">
            {text}
        </h{level}>
        """,
        unsafe_allow_html=True
    )


def add_focus_management(element_id: str):
    """Add focus management for dynamic content.
    
    Args:
        element_id: ID of element to focus
    """
    st.markdown(
        f"""
        <script>
            setTimeout(() => {{
                const element = document.getElementById('{element_id}');
                if (element) {{
                    element.focus();
                    element.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
                }}
            }}, 100);
        </script>
        """,
        unsafe_allow_html=True
    )


def render_status_with_icon(
    status: str, 
    message: str, 
    use_color: bool = True,
    icon_only: bool = False
):
    """Render status message with both color and icon for accessibility.
    
    Args:
        status: Status type ('success', 'error', 'warning', 'info')
        message: Status message
        use_color: Whether to use color styling
        icon_only: Whether to show only icon (for space-constrained areas)
    """
    status_config = {
        'success': {'icon': '✅', 'color': '#28a745', 'bg': '#d4edda'},
        'error': {'icon': '❌', 'color': '#dc3545', 'bg': '#f8d7da'},
        'warning': {'icon': '⚠️', 'color': '#ffc107', 'bg': '#fff3cd'},
        'info': {'icon': 'ℹ️', 'color': '#17a2b8', 'bg': '#d1ecf1'}
    }
    
    config = status_config.get(status, status_config['info'])
    
    if icon_only:
        st.markdown(
            f"""
            <span role="img" 
                  aria-label="{status}" 
                  title="{message}"
                  style="font-size: 1.2em;">
                {config['icon']}
            </span>
            """,
            unsafe_allow_html=True
        )
    else:
        color_style = f"color: {config['color']}; background-color: {config['bg']}; padding: 0.5rem; border-radius: 0.25rem;" if use_color else ""
        
        st.markdown(
            f"""
            <div role="status" 
                 aria-live="polite"
                 style="{color_style} margin: 0.5rem 0;">
                <span role="img" aria-label="{status}" style="margin-right: 0.5rem;">
                    {config['icon']}
                </span>
                {message}
            </div>
            """,
            unsafe_allow_html=True
        )


def add_table_accessibility(table_caption: str, table_id: Optional[str] = None):
    """Add accessibility features to tables.
    
    Args:
        table_caption: Caption describing the table content
        table_id: Optional ID for the table
    """
    table_id = table_id or f"table_{uuid.uuid4().hex[:8]}"
    
    st.markdown(
        f"""
        <div role="table" aria-label="{table_caption}" id="{table_id}">
            <div role="caption" style="font-weight: bold; margin-bottom: 0.5rem;">
                {table_caption}
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )


def ensure_minimum_contrast():
    """Add CSS to ensure minimum color contrast for WCAG AA compliance."""
    st.markdown(
        """
        <style>
            /* Ensure minimum contrast for text elements */
            .stMarkdown p, .stMarkdown li, .stMarkdown span {
                color: #212529 !important; /* High contrast dark text */
            }
            
            /* Ensure button contrast */
            .stButton > button {
                background-color: #0066cc !important;
                color: #ffffff !important;
                border: 2px solid #004499 !important;
            }
            
            .stButton > button:hover {
                background-color: #004499 !important;
                border-color: #002266 !important;
            }
            
            .stButton > button:focus {
                outline: 3px solid #ffbf00 !important;
                outline-offset: 2px !important;
            }
            
            /* Ensure selectbox contrast */
            .stSelectbox > div > div {
                background-color: #ffffff !important;
                color: #212529 !important;
                border: 2px solid #6c757d !important;
            }
            
            /* Focus indicators for all interactive elements */
            .stSelectbox > div > div:focus-within,
            .stTextInput > div > div > input:focus,
            .stTextArea > div > div > textarea:focus {
                outline: 3px solid #ffbf00 !important;
                outline-offset: 2px !important;
                border-color: #0066cc !important;
            }
            
            /* Ensure sidebar navigation contrast */
            .css-1d391kg {
                background-color: #f8f9fa !important;
            }
            
            .css-1d391kg .stMarkdown {
                color: #212529 !important;
            }
        </style>
        """,
        unsafe_allow_html=True
    )
