"""Alignment file management and visualization page."""

from __future__ import annotations

import os

import pandas as pd
import plotly.express as px
import requests
import streamlit as st

from scripts.dashboard.auth import require_auth


API_BASE = os.getenv("API_BASE_URL", "http://localhost:8000/api/v1")


def render_alignments_page() -> None:
    """Render the Alignment Files page."""
    require_auth()
    
    st.title("🧬 Alignment Files")
    st.markdown("Upload, browse, and visualize BAM/CRAM alignment files")
    
    tab_browse, tab_upload, tab_view, tab_stats = st.tabs([
        "📁 Browse", "📤 Upload", "🔬 View", "📊 Stats"
    ])
    
    with tab_browse:
        _render_browse_tab()
    
    with tab_upload:
        _render_upload_tab()
    
    with tab_view:
        _render_view_tab()
    
    with tab_stats:
        _render_stats_tab()


def _get_auth_headers() -> dict:
    """Get authentication headers from session."""
    token = st.session_state.get("access_token", "")
    return {"Authorization": f"Bearer {token}"}


def _render_browse_tab() -> None:
    """Render the Browse tab."""
    st.subheader("Browse Alignment Files")
    
    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        has_index_filter = st.selectbox(
            "Index Status",
            options=[None, True, False],
            format_func=lambda x: "All" if x is None else ("Has Index" if x else "No Index"),
            key="browse_index_filter",
        )
    with col2:
        limit = st.slider("Max Results", 10, 200, 50, key="browse_limit")
    with col3:
        if st.button("🔄 Refresh", key="browse_refresh"):
            st.rerun()
    
    # Fetch alignments
    try:
        params = {"limit": limit}
        if has_index_filter is not None:
            params["has_index"] = has_index_filter
        
        response = requests.get(
            f"{API_BASE}/genomics/alignments",
            headers=_get_auth_headers(),
            params=params,
            timeout=30,
        )
        
        if response.status_code == 200:
            alignments = response.json()
            
            if alignments:
                df = pd.DataFrame(alignments)
                df["created_at"] = pd.to_datetime(df["created_at"]).dt.strftime("%Y-%m-%d %H:%M")
                
                # Display table
                st.dataframe(
                    df[["filename", "file_format", "has_index", "reference_genome", "total_reads", "created_at"]],
                    use_container_width=True,
                    hide_index=True,
                )
                
                # Store selected alignment for other tabs
                selected = st.selectbox(
                    "Select alignment for viewing",
                    options=df["id"].tolist(),
                    format_func=lambda x: df[df["id"] == x]["filename"].values[0],
                    key="selected_alignment",
                )
                st.session_state["selected_alignment_id"] = selected
            else:
                st.info("No alignment files found. Upload one in the Upload tab.")
        else:
            st.error(f"Error fetching alignments: {response.status_code}")
    
    except Exception as e:
        st.error(f"Error: {e}")


def _render_upload_tab() -> None:
    """Render the Upload tab."""
    st.subheader("Upload Alignment File")
    
    st.markdown("""
    Upload BAM or CRAM alignment files. For region queries, an index file (.bai/.crai) is required.
    
    **Supported formats:**
    - BAM files (.bam) with optional .bai index
    - CRAM files (.cram) with optional .crai index
    
    **Size limit:** 500MB
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        alignment_file = st.file_uploader(
            "Alignment File (BAM/CRAM)",
            type=["bam", "cram"],
            key="upload_alignment",
        )
    
    with col2:
        index_file = st.file_uploader(
            "Index File (optional)",
            type=["bai", "crai"],
            key="upload_index",
        )
    
    # Optional experiment linking
    experiment_id = st.text_input(
        "Experiment ID (optional)",
        placeholder="UUID of experiment to link",
        key="upload_experiment_id",
    )
    
    if st.button("📤 Upload", key="upload_submit", type="primary"):
        if not alignment_file:
            st.error("Please select an alignment file")
            return
        
        with st.spinner("Uploading and processing..."):
            try:
                files = {"file": (alignment_file.name, alignment_file.getvalue())}
                if index_file:
                    files["index_file"] = (index_file.name, index_file.getvalue())
                
                params = {}
                if experiment_id:
                    params["experiment_id"] = experiment_id
                
                response = requests.post(
                    f"{API_BASE}/genomics/alignments/upload",
                    headers=_get_auth_headers(),
                    files=files,
                    params=params,
                    timeout=120,
                )
                
                if response.status_code == 200:
                    result = response.json()
                    if result.get("success"):
                        st.success(f"✅ Uploaded {result['filename']} successfully!")
                        st.json(result.get("stats", {}))
                        st.balloons()
                    else:
                        st.error(f"Upload failed: {result.get('error', 'Unknown error')}")
                else:
                    st.error(f"Upload failed: {response.status_code} - {response.text}")
            
            except Exception as e:
                st.error(f"Upload error: {e}")


def _render_view_tab() -> None:
    """Render the View tab for region queries."""
    st.subheader("View Alignment Reads")
    
    alignment_id = st.session_state.get("selected_alignment_id")
    
    if not alignment_id:
        st.info("Select an alignment file in the Browse tab first.")
        return
    
    # Region input
    col1, col2 = st.columns([3, 1])
    
    with col1:
        region = st.text_input(
            "Genomic Region",
            placeholder="chr1:1000-2000",
            help="Format: chromosome:start-end (e.g., chr1:1000000-1001000)",
            key="view_region",
        )
    
    with col2:
        page_size = st.selectbox(
            "Reads per page",
            options=[50, 100, 200, 500],
            index=1,
            key="view_page_size",
        )
    
    col3, col4 = st.columns(2)
    
    with col3:
        if st.button("🔍 Fetch Reads", key="view_fetch", type="primary"):
            if not region:
                st.error("Please enter a genomic region")
                return
            
            with st.spinner("Fetching reads..."):
                try:
                    response = requests.get(
                        f"{API_BASE}/genomics/alignments/{alignment_id}/reads",
                        headers=_get_auth_headers(),
                        params={"region": region, "page_size": page_size},
                        timeout=60,
                    )
                    
                    if response.status_code == 200:
                        result = response.json()
                        if result.get("success"):
                            reads = result.get("reads", [])
                            st.success(f"Found {len(reads)} reads in region {region}")
                            
                            if reads:
                                df = pd.DataFrame(reads)
                                st.dataframe(
                                    df[["query_name", "reference_start", "mapping_quality", "cigar", "is_reverse"]],
                                    use_container_width=True,
                                    hide_index=True,
                                )
                        else:
                            st.error(result.get("error", "Failed to fetch reads"))
                    elif response.status_code == 400:
                        st.error("Index file required for region queries. Upload an index file.")
                    else:
                        st.error(f"Error: {response.status_code}")
                
                except Exception as e:
                    st.error(f"Error: {e}")
    
    with col4:
        if st.button("📊 Show Coverage", key="view_coverage"):
            if not region:
                st.error("Please enter a genomic region")
                return
            
            with st.spinner("Calculating coverage..."):
                try:
                    response = requests.get(
                        f"{API_BASE}/genomics/alignments/{alignment_id}/coverage",
                        headers=_get_auth_headers(),
                        params={"region": region, "bin_size": 100},
                        timeout=60,
                    )
                    
                    if response.status_code == 200:
                        result = response.json()
                        if result.get("success"):
                            coverage = result.get("coverage", [])
                            
                            if coverage:
                                df = pd.DataFrame(coverage)
                                fig = px.line(
                                    df, x="position", y="coverage",
                                    title=f"Coverage: {region}",
                                    labels={"position": "Position", "coverage": "Coverage Depth"},
                                )
                                st.plotly_chart(fig, use_container_width=True)
                        else:
                            st.error(result.get("error", "Failed to get coverage"))
                    else:
                        st.error(f"Error: {response.status_code}")
                
                except Exception as e:
                    st.error(f"Error: {e}")


def _render_stats_tab() -> None:
    """Render the Stats tab."""
    st.subheader("Alignment Statistics")
    
    alignment_id = st.session_state.get("selected_alignment_id")
    
    if not alignment_id:
        st.info("Select an alignment file in the Browse tab first.")
        return
    
    try:
        response = requests.get(
            f"{API_BASE}/genomics/alignments/{alignment_id}",
            headers=_get_auth_headers(),
            timeout=30,
        )
        
        if response.status_code == 200:
            data = response.json()
            
            # File info
            st.markdown("### File Information")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Filename", data.get("filename", "Unknown"))
            with col2:
                st.metric("Format", data.get("file_format", "Unknown"))
            with col3:
                st.metric("Has Index", "✅ Yes" if data.get("has_index") else "❌ No")
            
            # Reference info
            st.markdown("### Reference Genome")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Reference", data.get("reference_genome") or "Unknown")
            with col2:
                st.metric("Chromosomes", data.get("num_references") or 0)
            
            # Read statistics
            st.markdown("### Read Statistics")
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                total = data.get("total_reads") or 0
                st.metric("Total Reads", f"{total:,}")
            with col2:
                mapped = data.get("mapped_reads") or 0
                st.metric("Mapped Reads", f"{mapped:,}")
            with col3:
                unmapped = data.get("unmapped_reads") or 0
                st.metric("Unmapped Reads", f"{unmapped:,}")
            with col4:
                dup_rate = data.get("duplicate_rate") or 0
                st.metric("Duplicate Rate", f"{dup_rate:.2%}")
            
            # Read groups
            read_groups = data.get("read_groups")
            if read_groups:
                st.markdown("### Read Groups")
                st.dataframe(pd.DataFrame(read_groups), use_container_width=True, hide_index=True)
        
        else:
            st.error(f"Error fetching alignment details: {response.status_code}")
    
    except Exception as e:
        st.error(f"Error: {e}")
