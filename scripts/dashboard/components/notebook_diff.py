"""Notebook diff viewer component for review comparisons."""

from typing import Dict, List

import httpx
import streamlit as st


def render_notebook_diff(review_id: str, api_base: str = "http://localhost:8000") -> None:
    """Render the notebook diff view for a review.
    
    Args:
        review_id: UUID of the notebook review
        api_base: Base URL for the API server
    """
    st.subheader("📋 Notebook Changes")
    
    # Load diff data
    try:
        with st.spinner("Computing notebook diff..."):
            response = httpx.get(f"{api_base}/api/v1/reviews/{review_id}/diff")
            response.raise_for_status()
            diff_data = response.json()
    except httpx.ConnectError:
        st.warning("⚠️ API server unavailable. Showing demo data.")
        diff_data = _get_demo_diff()
    except Exception as e:
        st.error(f"Failed to load diff: {e}")
        diff_data = {"added": [], "removed": [], "modified": [], "unchanged": []}
    
    # Extract change categories
    added = diff_data.get("added", [])
    removed = diff_data.get("removed", [])
    modified = diff_data.get("modified", [])
    unchanged = diff_data.get("unchanged", [])
    
    # Summary statistics
    total_changes = len(added) + len(removed) + len(modified)
    
    if total_changes == 0:
        st.info("✨ No changes detected. The notebook matches the reviewed snapshot.")
        return
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Added", len(added), delta=len(added) if added else None)
    with col2:
        st.metric("Removed", len(removed), delta=-len(removed) if removed else None)
    with col3:
        st.metric("Modified", len(modified), delta=len(modified) if modified else None)
    with col4:
        st.metric("Unchanged", len(unchanged))
    
    st.divider()
    
    # Change sections
    if added:
        _render_added_cells(added)
    
    if removed:
        _render_removed_cells(removed)
    
    if modified:
        _render_modified_cells(modified)
    
    if unchanged:
        _render_unchanged_cells(unchanged)


def _render_added_cells(added: List[Dict]) -> None:
    """Render added cells section."""
    with st.expander(f"➕ **Added Cells ({len(added)})**", expanded=True):
        st.success(f"Found {len(added)} new cell(s)")
        
        for cell in added:
            cell_type = cell.get("cell_type", "unknown")
            cell_index = cell.get("index", "?")
            source_preview = cell.get("source_preview", "")
            
            # Cell header
            cell_badge = _get_cell_type_badge(cell_type)
            st.markdown(f"**Cell #{cell_index}** {cell_badge}")
            
            # Cell content
            if source_preview:
                language = "python" if cell_type == "code" else "markdown"
                st.code(source_preview, language=language)
            else:
                st.caption("_Empty cell_")
            
            st.divider()


def _render_removed_cells(removed: List[Dict]) -> None:
    """Render removed cells section."""
    with st.expander(f"➖ **Removed Cells ({len(removed)})**", expanded=True):
        st.error(f"Found {len(removed)} deleted cell(s)")
        
        for cell in removed:
            cell_type = cell.get("cell_type", "unknown")
            cell_index = cell.get("index", "?")
            source_preview = cell.get("source_preview", "")
            
            # Cell header
            cell_badge = _get_cell_type_badge(cell_type)
            st.markdown(f"**Cell #{cell_index}** {cell_badge}")
            
            # Cell content
            if source_preview:
                language = "python" if cell_type == "code" else "markdown"
                st.code(source_preview, language=language)
            else:
                st.caption("_Empty cell_")
            
            st.divider()


def _render_modified_cells(modified: List[Dict]) -> None:
    """Render modified cells section with side-by-side comparison."""
    with st.expander(f"🔄 **Modified Cells ({len(modified)})**", expanded=True):
        st.warning(f"Found {len(modified)} changed cell(s)")
        
        for cell in modified:
            cell_type = cell.get("cell_type", "unknown")
            cell_index = cell.get("index", "?")
            old_preview = cell.get("old_source_preview", "")
            new_preview = cell.get("new_source_preview", "")
            
            # Cell header
            cell_badge = _get_cell_type_badge(cell_type)
            st.markdown(f"**Cell #{cell_index}** {cell_badge}")
            
            # Side-by-side comparison
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Before:**")
                if old_preview:
                    language = "python" if cell_type == "code" else "markdown"
                    st.code(old_preview, language=language)
                else:
                    st.caption("_Empty cell_")
            
            with col2:
                st.markdown("**After:**")
                if new_preview:
                    language = "python" if cell_type == "code" else "markdown"
                    st.code(new_preview, language=language)
                else:
                    st.caption("_Empty cell_")
            
            st.divider()


def _render_unchanged_cells(unchanged: List[Dict]) -> None:
    """Render unchanged cells section (collapsed by default)."""
    with st.expander(f"✅ **Unchanged Cells ({len(unchanged)})**", expanded=False):
        st.info(f"Found {len(unchanged)} unchanged cell(s)")
        
        # Group by cell type for summary
        cell_types = {}
        for cell in unchanged:
            cell_type = cell.get("cell_type", "unknown")
            cell_types[cell_type] = cell_types.get(cell_type, 0) + 1
        
        # Summary by type
        for cell_type, count in cell_types.items():
            badge = _get_cell_type_badge(cell_type)
            st.markdown(f"- {count} {badge} cells")
        
        # Option to show all unchanged cells
        if st.checkbox("Show all unchanged cells", key="show_unchanged"):
            for cell in unchanged:
                cell_type = cell.get("cell_type", "unknown")
                cell_index = cell.get("index", "?")
                source_preview = cell.get("source_preview", "")
                
                # Cell header
                cell_badge = _get_cell_type_badge(cell_type)
                st.markdown(f"**Cell #{cell_index}** {cell_badge}")
                
                # Cell content (truncated for unchanged)
                if source_preview:
                    # Truncate long content for unchanged cells
                    preview = source_preview[:200] + "..." if len(source_preview) > 200 else source_preview
                    language = "python" if cell_type == "code" else "markdown"
                    st.code(preview, language=language)
                else:
                    st.caption("_Empty cell_")


def _get_cell_type_badge(cell_type: str) -> str:
    """Get a colored badge for cell type."""
    badges = {
        "code": ":blue[Code]",
        "markdown": ":green[Markdown]",
        "raw": ":gray[Raw]",
    }
    return badges.get(cell_type, f":orange[{cell_type.title()}]")


def _get_demo_diff() -> Dict:
    """Return demo diff data when API is unavailable."""
    return {
        "added": [
            {
                "index": 0,
                "cell_type": "markdown",
                "source_preview": "# New Analysis Section\n\nThis section contains additional analysis that was added after the initial review."
            },
            {
                "index": 5,
                "cell_type": "code",
                "source_preview": "# New visualization\nimport matplotlib.pyplot as plt\n\nplt.figure(figsize=(10, 6))\nplt.plot(data['x'], data['y'])\nplt.title('New Plot')\nplt.show()"
            }
        ],
        "removed": [
            {
                "index": 3,
                "cell_type": "code",
                "source_preview": "# Deprecated code that was removed\nold_function(deprecated_params)"
            }
        ],
        "modified": [
            {
                "index": 2,
                "cell_type": "code",
                "old_source_preview": "# Original implementation\nresult = simple_calculation(x, y)",
                "new_source_preview": "# Improved implementation with error handling\ntry:\n    result = enhanced_calculation(x, y, validate=True)\nexcept ValueError as e:\n    print(f'Calculation error: {e}')\n    result = None"
            },
            {
                "index": 4,
                "cell_type": "markdown",
                "old_source_preview": "## Results\n\nThe analysis shows basic trends.",
                "new_source_preview": "## Results\n\nThe analysis shows significant trends with statistical confidence (p < 0.05).\n\n### Key Findings:\n- Finding 1\n- Finding 2"
            }
        ],
        "unchanged": [
            {
                "index": 1,
                "cell_type": "code",
                "source_preview": "import pandas as pd\nimport numpy as np\n\n# Load data\ndata = pd.read_csv('dataset.csv')"
            },
            {
                "index": 6,
                "cell_type": "markdown",
                "source_preview": "## Conclusion\n\nThis notebook demonstrates the analysis workflow."
            }
        ]
    }
