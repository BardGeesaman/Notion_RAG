"""Theme management for Streamlit dashboard."""
from __future__ import annotations

from typing import List

import streamlit as st

THEMES = {
    "dark": {
        "background": "#0E1117",
        "text": "#FAFAFA",
        "sidebar_background": "#262730",
        "sidebar_text": "#FAFAFA",
        "primary": "#FF4B4B",
        "secondary": "#262730",
        "accent": "#00D4AA",
        "success": "#00D4AA",
        "warning": "#FFA726",
        "error": "#FF4B4B",
        "info": "#1E88E5",
    },
    "light": {
        "background": "#FFFFFF",
        "text": "#262730",
        "sidebar_background": "#F0F2F6",
        "sidebar_text": "#262730",
        "primary": "#FF4B4B",
        "secondary": "#F0F2F6",
        "accent": "#00D4AA",
        "success": "#00D4AA",
        "warning": "#FFA726",
        "error": "#FF4B4B",
        "info": "#1E88E5",
    },
}


def apply_theme(theme_name: str) -> None:
    """
    Apply a theme by injecting CSS into the Streamlit app.

    Args:
        theme_name: Name of the theme ("dark" or "light")
    """
    if theme_name not in THEMES:
        theme_name = "light"  # Default to light

    theme = THEMES[theme_name]

    css = f"""
    <style>
    :root {{
        --background-color: {theme['background']};
        --text-color: {theme['text']};
        --sidebar-background: {theme['sidebar_background']};
        --sidebar-text: {theme['sidebar_text']};
        --primary-color: {theme['primary']};
        --secondary-color: {theme['secondary']};
        --accent-color: {theme['accent']};
    }}

    /* Override Streamlit default styles */
    .stApp {{
        background-color: {theme['background']};
        color: {theme['text']};
    }}

    /* Main content area */
    .main .block-container {{
        background-color: {theme['background']};
        color: {theme['text']};
    }}

    /* Sidebar */
    [data-testid="stSidebar"] {{
        background-color: {theme['sidebar_background']} !important;
    }}

    [data-testid="stSidebar"] * {{
        color: {theme['sidebar_text']} !important;
    }}

    /* Text elements */
    h1, h2, h3, h4, h5, h6, p, div, span, label {{
        color: {theme['text']} !important;
    }}

    /* Buttons */
    .stButton > button {{
        background-color: {theme['primary']};
        color: white;
        border-radius: 0.25rem;
    }}

    .stButton > button:hover {{
        background-color: {theme['accent']};
    }}

    /* Input fields */
    .stTextInput > div > div > input {{
        background-color: {theme['background']};
        color: {theme['text']};
        border-color: {theme['secondary']};
    }}

    .stSelectbox > div > div > select {{
        background-color: {theme['background']};
        color: {theme['text']};
    }}

    /* Dataframes */
    .stDataFrame {{
        background-color: {theme['background']};
        color: {theme['text']};
    }}

    /* Expanders */
    .streamlit-expanderHeader {{
        background-color: {theme['secondary']};
        color: {theme['text']};
    }}

    /* Code blocks */
    .stCodeBlock {{
        background-color: {theme['secondary']};
        color: {theme['text']};
    }}
    </style>
    """

    st.markdown(css, unsafe_allow_html=True)


def get_theme_names() -> List[str]:
    """
    Get list of available theme names.

    Returns:
        List of theme names
    """
    return list(THEMES.keys())
