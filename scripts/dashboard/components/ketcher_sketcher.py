"""
Ketcher chemical structure editor component for Streamlit.

Embeds Ketcher (open-source chemical editor) via iframe with
bidirectional SMILES communication using window.postMessage.
"""

from __future__ import annotations

import uuid
from typing import Optional

import streamlit as st
import streamlit.components.v1 as components


def build_ketcher_html(
    initial_smiles: Optional[str] = None,
    width: int = 800,
    height: int = 600,
    component_key: Optional[str] = None,
) -> str:
    """
    Build standalone HTML for Ketcher editor with SMILES I/O.

    Args:
        initial_smiles: SMILES string to load initially
        width: Component width in pixels
        height: Component height in pixels
        component_key: Unique key for this component instance

    Returns:
        HTML string with embedded Ketcher iframe
    """
    dom_id = f"ketcher_{component_key or uuid.uuid4().hex}"
    initial_smiles_escaped = (initial_smiles or "").replace("'", "\\'").replace('"', '\\"')
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <style>
        body {{ margin: 0; padding: 0; overflow: hidden; }}
        #ketcher-frame {{ width: 100%; height: {height}px; border: 1px solid #ccc; }}
        #controls {{ padding: 10px; background: #f5f5f5; }}
        button {{ margin: 5px; padding: 8px 16px; cursor: pointer; }}
    </style>
</head>
<body>
    <div id="controls">
        <button onclick="exportSMILES()">Get SMILES</button>
        <span id="smiles-output" style="margin-left: 20px; font-family: monospace;"></span>
    </div>
    <iframe id="ketcher-frame" src="https://unpkg.com/ketcher-standalone@2.7.2/dist/index.html"></iframe>
    
    <script>
        let ketcherFrame = document.getElementById('ketcher-frame');
        let ketcherInstance = null;
        
        // Wait for Ketcher to load
        ketcherFrame.onload = function() {{
            ketcherInstance = ketcherFrame.contentWindow.ketcher;
            
            // Load initial SMILES if provided
            const initialSmiles = '{initial_smiles_escaped}';
            if (initialSmiles && ketcherInstance) {{
                setTimeout(() => {{
                    ketcherInstance.setMolecule(initialSmiles);
                }}, 500);
            }}
        }};
        
        // Export SMILES function
        async function exportSMILES() {{
            if (!ketcherInstance) {{
                alert('Ketcher not loaded yet');
                return;
            }}
            
            try {{
                const smiles = await ketcherInstance.getSmiles();
                document.getElementById('smiles-output').textContent = smiles || '(empty)';
                
                // Send SMILES to Streamlit via parent window
                if (window.parent) {{
                    window.parent.postMessage({{
                        type: 'ketcher-smiles',
                        key: '{dom_id}',
                        smiles: smiles
                    }}, '*');
                }}
            }} catch (err) {{
                console.error('Failed to export SMILES:', err);
                document.getElementById('smiles-output').textContent = 'Error: ' + err.message;
            }}
        }}
    </script>
</body>
</html>
"""
    return html


def render_ketcher_editor(
    key: str = "ketcher",
    initial_smiles: Optional[str] = None,
    width: int = 800,
    height: int = 600,
) -> None:
    """
    Render Ketcher chemical structure editor in Streamlit.

    Args:
        key: Unique key for this component instance
        initial_smiles: SMILES string to load initially
        width: Component width in pixels
        height: Component height in pixels

    Usage:
        render_ketcher_editor(key="my_ketcher", initial_smiles="CCO")
        smiles = get_smiles_from_session("my_ketcher")
    """
    html = build_ketcher_html(
        initial_smiles=initial_smiles,
        width=width,
        height=height,
        component_key=key,
    )
    
    # Render component
    components.html(html, height=height + 60, scrolling=False)
    
    # Note: SMILES extraction via postMessage requires JavaScript bridge
    # For now, users click "Get SMILES" button which displays in component
    # Full Streamlit integration would require custom component package


def get_smiles_from_session(key: str) -> Optional[str]:
    """
    Get SMILES from Streamlit session state.

    Args:
        key: Component key used in render_ketcher_editor

    Returns:
        SMILES string if available, None otherwise

    Note: Due to iframe isolation, SMILES must be manually copied or
    requires custom Streamlit component for automatic extraction.
    For MVP, users can copy SMILES from component display.
    """
    session_key = f"ketcher_smiles_{key}"
    return st.session_state.get(session_key)

