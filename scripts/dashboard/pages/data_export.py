"""Data Export Wizard page."""

from __future__ import annotations

import os

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_data_export_page() -> None:
    """Render the Data Export Wizard page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📦 Data Export Wizard")
    st.caption("Export datasets, experiments, and compounds to CSV, Excel, or JSON.")
    
    tab1, tab2, tab3 = st.tabs(["Quick Export", "Package Builder", "Export History"])
    
    with tab1:
        render_quick_export_tab()
    
    with tab2:
        render_package_builder_tab()
    
    with tab3:
        render_export_history_tab()


def render_quick_export_tab() -> None:
    """Render quick export tab."""
    st.subheader("Quick Export")
    
    entity_type = st.selectbox(
        "Entity Type",
        options=["Dataset", "Experiment", "Compound"],
        index=0,
    )
    
    # Entity ID input (simplified - would be dropdown in production)
    entity_id = st.text_input(
        f"{entity_type} ID",
        placeholder=f"UUID of {entity_type.lower()}",
    )
    
    format_option = st.selectbox(
        "Format",
        options=["CSV", "Excel", "JSON"],
        index=0,
    )
    
    if entity_type == "Dataset":
        include_metadata = st.checkbox("Include metadata header", value=True)
    
    if st.button("Download", type="primary", disabled=not entity_id):
        try:
            if entity_type == "Dataset":
                url = f"{API_BASE}/api/v1/export/dataset/{entity_id}?format={format_option.lower()}&include_metadata={include_metadata}"
            elif entity_type == "Experiment":
                url = f"{API_BASE}/api/v1/export/experiment/{entity_id}?format={format_option.lower()}"
            else:
                st.warning("Compound export uses package builder")
                return
            
            with httpx.Client(timeout=60) as client:
                response = client.get(url)
                response.raise_for_status()
                
                # Offer download
                filename = f"{entity_type.lower()}_{entity_id[:8]}.{format_option.lower()}"
                st.download_button(
                    "💾 Save File",
                    data=response.content,
                    file_name=filename,
                    mime=response.headers.get("content-type", "application/octet-stream"),
                )
                
                # Track in history
                if "export_history" not in st.session_state:
                    st.session_state["export_history"] = []
                
                st.session_state["export_history"].append({
                    "entity": f"{entity_type} {entity_id[:8]}",
                    "format": format_option,
                    "size": len(response.content),
                })
                
        except Exception as e:
            st.error(f"Export failed: {e}")


def render_package_builder_tab() -> None:
    """Render package builder tab."""
    st.subheader("Package Builder")
    st.caption("Bundle multiple exports into a ZIP package")
    
    # Package list
    if "package_items" not in st.session_state:
        st.session_state["package_items"] = []
    
    # Add item
    col1, col2, col3 = st.columns([2, 1, 1])
    
    with col1:
        item_id = st.text_input("Item ID", placeholder="UUID to add to package", key="package_item_id")
    
    with col2:
        if st.button("➕ Add to Package"):
            if item_id and item_id not in st.session_state["package_items"]:
                st.session_state["package_items"].append(item_id)
    
    with col3:
        if st.button("🗑️ Clear All"):
            st.session_state["package_items"] = []
    
    # Show current package
    if st.session_state["package_items"]:
        st.markdown(f"**Package Contents** ({len(st.session_state['package_items'])} items):")
        for item in st.session_state["package_items"]:
            st.code(item, language=None)
        
        if st.button("📦 Generate Package", type="primary"):
            try:
                with httpx.Client(timeout=120) as client:
                    response = client.post(
                        f"{API_BASE}/api/v1/export/package",
                        json={
                            "items": st.session_state["package_items"],
                            "format": "zip",
                        },
                    )
                    response.raise_for_status()
                    
                    st.download_button(
                        "💾 Download ZIP Package",
                        data=response.content,
                        file_name="export_package.zip",
                        mime="application/zip",
                    )
                    
            except Exception as e:
                st.error(f"Package generation failed: {e}")
    else:
        st.info("Add items to package using the input above")


def render_export_history_tab() -> None:
    """Render export history tab."""
    st.subheader("Export History")
    
    history = st.session_state.get("export_history", [])
    
    if history:
        import pandas as pd
        
        df = pd.DataFrame(history)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        st.caption(f"Total exports: {len(history)}")
        
        if st.button("Clear History"):
            st.session_state["export_history"] = []
            st.rerun()
    else:
        st.info("No export history yet. Export some data from the Quick Export tab!")


if __name__ == "__main__":
    render_data_export_page()

