"""Imaging Browser dashboard page."""

from __future__ import annotations

import os
import tempfile
import zipfile
from typing import Any, Dict, List, Optional
import json

import httpx
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import numpy as np


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_imaging_browser_page() -> None:
    """Render the Imaging Browser dashboard."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🔬 Microscopy Image Browser")
    st.caption("Import, browse, and analyze high-content imaging data from multiple vendor platforms.")

    # Create tabs
    tab_import, tab_browse, tab_qc, tab_instruments = st.tabs([
        "📤 Batch Import", "🖼️ 5D Browser", "📊 QC Dashboard", "🔧 Instruments"
    ])

    with tab_import:
        render_batch_import_tab()
    
    with tab_browse:
        render_5d_browser_tab()
    
    with tab_qc:
        render_qc_dashboard_tab()
    
    with tab_instruments:
        render_instruments_tab()


def render_batch_import_tab() -> None:
    """Render the Batch Import tab."""
    st.subheader("Batch Import from Vendor Exports")
    st.caption("Import images from Opera, ImageXpress, and Cell Voyager export directories")
    
    # Import section
    st.markdown("### 📁 Import Configuration")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Import path input
        import_path = st.text_input(
            "Export Directory Path",
            placeholder="/path/to/vendor/export",
            help="Path to vendor export directory containing images and metadata"
        )
        
        # Vendor selection
        vendor = st.selectbox(
            "Vendor Format",
            ["Auto-detect", "Opera/Operetta", "ImageXpress", "Cell Voyager"],
            help="Select vendor format or use auto-detection"
        )
        
        # Plate association
        plate_id = st.text_input(
            "Plate ID (optional)",
            placeholder="Enter existing plate ID or leave blank",
            help="Associate with existing plate or create new one"
        )
        
        create_wells = st.checkbox(
            "Create missing wells",
            value=True,
            help="Automatically create wells that don't exist in the database"
        )
    
    with col2:
        st.markdown("### 📋 Supported Formats")
        st.info("""
        **Opera/Operetta:**
        - Index.xml + Images/
        - r{row}c{col}f{field}... naming
        
        **ImageXpress:**
        - HTD files + TimePoint directories
        - {well}_{site}_w{wavelength}.TIF
        
        **Cell Voyager:**
        - MLF files + .flex images
        - W{well}F{field}T{time}Z{z}C{channel}
        """)
    
    # Import button and status
    if st.button("🚀 Start Import", type="primary", disabled=not import_path):
        if import_path:
            start_batch_import(import_path, vendor, plate_id, create_wells)
    
    # Display recent imports
    st.markdown("### 📈 Recent Imports")
    display_import_history()


def render_5d_browser_tab() -> None:
    """Render the 5D Browser tab."""
    st.subheader("5D Image Browser")
    st.caption("Browse images across Plate, Well, Channel, Z-slice, and Timepoint dimensions")
    
    # Filters section
    st.markdown("### 🔍 Filters")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        # Plate selector
        plates = get_available_plates()
        selected_plate = st.selectbox(
            "Plate",
            options=[None] + plates,
            format_func=lambda x: "All plates" if x is None else f"Plate {x}",
            help="Select a specific plate or browse all"
        )
        
        # Channel selector
        channels = get_available_channels(selected_plate)
        selected_channels = st.multiselect(
            "Channels",
            options=channels,
            default=channels[:3] if len(channels) > 3 else channels,
            help="Select channels to display"
        )
    
    with col2:
        # Well position filter
        well_position = st.text_input(
            "Well Position",
            placeholder="e.g., A01, B12",
            help="Filter by specific well position"
        )
        
        # Z-slice range
        z_range = get_z_range(selected_plate)
        if z_range[1] > z_range[0]:
            z_slice = st.slider(
                "Z-slice",
                min_value=z_range[0],
                max_value=z_range[1],
                value=z_range[0],
                help="Select Z-slice to view"
            )
        else:
            z_slice = 0
    
    with col3:
        # Timepoint range
        t_range = get_t_range(selected_plate)
        if t_range[1] > t_range[0]:
            timepoint = st.slider(
                "Timepoint",
                min_value=t_range[0],
                max_value=t_range[1],
                value=t_range[0],
                help="Select timepoint to view"
            )
        else:
            timepoint = 0
        
        # QC filter
        qc_filter = st.selectbox(
            "QC Status",
            ["All images", "Passed QC only", "Failed QC only"],
            help="Filter by quality control status"
        )
    
    # Image grid
    st.markdown("### 🖼️ Image Grid")
    
    if selected_channels:
        images = browse_images(
            plate_id=selected_plate,
            well_position=well_position,
            channels=selected_channels,
            z_slice=z_slice,
            timepoint=timepoint,
            qc_filter=qc_filter
        )
        
        if images:
            display_image_grid(images)
        else:
            st.info("No images found matching the current filters")
    else:
        st.warning("Please select at least one channel to browse images")


def render_qc_dashboard_tab() -> None:
    """Render the QC Dashboard tab."""
    st.subheader("Quality Control Dashboard")
    st.caption("Monitor image quality metrics across plates and identify issues")
    
    # Plate selector
    plates = get_available_plates()
    selected_plate = st.selectbox(
        "Select Plate for QC Analysis",
        options=plates,
        help="Choose a plate to analyze quality metrics"
    )
    
    if selected_plate:
        # Get plate QC data
        qc_data = get_plate_qc_data(selected_plate)
        
        if qc_data:
            # Summary metrics
            st.markdown("### 📊 Quality Summary")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(
                    "Total Images",
                    qc_data.get("total_images", 0)
                )
            
            with col2:
                passed_count = qc_data.get("passed_count", 0)
                total_count = qc_data.get("total_images", 1)
                pass_rate = (passed_count / total_count) * 100
                st.metric(
                    "QC Pass Rate",
                    f"{pass_rate:.1f}%",
                    delta=f"{passed_count}/{total_count}"
                )
            
            with col3:
                avg_focus = qc_data.get("average_focus_score", 0)
                st.metric(
                    "Avg Focus Score",
                    f"{avg_focus:.3f}",
                    delta="Good" if avg_focus > 0.3 else "Poor"
                )
            
            with col4:
                failed_count = qc_data.get("failed_count", 0)
                st.metric(
                    "Failed Images",
                    failed_count,
                    delta="Issues detected" if failed_count > 0 else "All good"
                )
            
            # Focus heatmap
            st.markdown("### 🔥 Focus Quality Heatmap")
            focus_heatmap = create_focus_heatmap(qc_data.get("focus_heatmap", {}))
            if focus_heatmap:
                st.plotly_chart(focus_heatmap, use_container_width=True)
            
            # Issues summary
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("### ⚠️ Saturation Alerts")
                saturation_alerts = qc_data.get("saturation_alerts", [])
                if saturation_alerts:
                    for alert in saturation_alerts[:10]:  # Show top 10
                        st.error(f"Well {alert}: High saturation detected")
                else:
                    st.success("No saturation issues detected")
            
            with col2:
                st.markdown("### 💡 Uniformity Issues")
                uniformity_issues = qc_data.get("uniformity_issues", [])
                if uniformity_issues:
                    for issue in uniformity_issues[:10]:  # Show top 10
                        st.warning(f"Well {issue}: Illumination non-uniformity")
                else:
                    st.success("No uniformity issues detected")
            
            # Recommendations
            st.markdown("### 💡 Recommendations")
            recommendations = qc_data.get("recommendations", [])
            if recommendations:
                for rec in recommendations:
                    st.info(f"💡 {rec}")
            else:
                st.success("No specific recommendations - image quality looks good!")
        
        else:
            st.warning("No QC data available for this plate")
    else:
        st.info("Please select a plate to view QC metrics")


def render_instruments_tab() -> None:
    """Render the Instruments tab."""
    st.subheader("Microscope Instrument Registry")
    st.caption("Manage microscopes, objectives, and channel configurations")
    
    # Instruments list
    st.markdown("### 🔬 Registered Microscopes")
    
    instruments = get_instruments()
    
    if instruments:
        # Create instruments table
        instruments_df = pd.DataFrame([
            {
                "Name": inst["name"],
                "Manufacturer": inst["manufacturer"],
                "Model": inst["model"],
                "Serial": inst.get("serial_number", "N/A"),
                "Location": inst.get("facility_location", "N/A"),
                "Objectives": len(inst.get("objectives", [])),
                "Channels": len(inst.get("channels", [])),
                "Status": "Active" if inst.get("is_active", True) else "Inactive"
            }
            for inst in instruments
        ])
        
        st.dataframe(instruments_df, use_container_width=True)
        
        # Instrument details
        st.markdown("### 🔍 Instrument Details")
        
        selected_instrument = st.selectbox(
            "Select Instrument",
            options=range(len(instruments)),
            format_func=lambda i: instruments[i]["name"]
        )
        
        if selected_instrument is not None:
            instrument = instruments[selected_instrument]
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("#### 🎯 Objectives")
                objectives = instrument.get("objectives", [])
                if objectives:
                    obj_df = pd.DataFrame([
                        {
                            "Name": obj["name"],
                            "Magnification": f"{obj['magnification']}x",
                            "NA": obj["numerical_aperture"],
                            "Immersion": obj["immersion"]
                        }
                        for obj in objectives
                    ])
                    st.dataframe(obj_df, use_container_width=True)
                else:
                    st.info("No objectives configured")
            
            with col2:
                st.markdown("#### 🌈 Channels")
                channels = instrument.get("channels", [])
                if channels:
                    ch_df = pd.DataFrame([
                        {
                            "Channel": ch["channel_name"],
                            "Fluorophore": ch.get("fluorophore", "N/A"),
                            "Exposure (ms)": ch.get("default_exposure_ms", "N/A"),
                            "Gain": ch.get("default_gain", "N/A")
                        }
                        for ch in channels
                    ])
                    st.dataframe(ch_df, use_container_width=True)
                else:
                    st.info("No channels configured")
    
    else:
        st.info("No instruments registered")
    
    # Add new instrument
    st.markdown("### ➕ Add New Microscope")
    
    with st.expander("Register New Instrument"):
        with st.form("new_instrument_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                name = st.text_input("Instrument Name", placeholder="e.g., Nikon Ti2-E #1")
                manufacturer = st.text_input("Manufacturer", placeholder="e.g., Nikon")
                model = st.text_input("Model", placeholder="e.g., Ti2-E")
            
            with col2:
                serial_number = st.text_input("Serial Number", placeholder="Optional")
                facility_location = st.text_input("Location", placeholder="e.g., Lab A, Room 101")
            
            submitted = st.form_submit_button("Register Instrument")
            
            if submitted and name and manufacturer and model:
                register_instrument(name, manufacturer, model, serial_number, facility_location)


# Helper functions for API calls and data processing

def start_batch_import(import_path: str, vendor: str, plate_id: Optional[str], create_wells: bool) -> None:
    """Start a batch import job."""
    try:
        # Map vendor selection
        vendor_map = {
            "Auto-detect": None,
            "Opera/Operetta": "opera",
            "ImageXpress": "imagexpress",
            "Cell Voyager": "cellvoyager"
        }
        
        data = {
            "import_path": import_path,
            "vendor": vendor_map.get(vendor),
            "create_missing_wells": create_wells
        }
        
        if plate_id:
            data["plate_id"] = plate_id
        
        response = httpx.post(f"{API_BASE}/imaging/import/batch", data=data)
        
        if response.status_code == 202:
            result = response.json()
            st.success(f"Import started! Job ID: {result['fileset_id']}")
            st.info(f"Vendor: {result['vendor']}, Total images: {result['total_images']}")
            
            # Store job ID in session state for tracking
            if "import_jobs" not in st.session_state:
                st.session_state.import_jobs = []
            st.session_state.import_jobs.append(result['fileset_id'])
            
        else:
            st.error(f"Import failed: {response.text}")
            
    except Exception as e:
        st.error(f"Error starting import: {str(e)}")


def display_import_history() -> None:
    """Display recent import jobs and their status."""
    if "import_jobs" in st.session_state and st.session_state.import_jobs:
        for job_id in st.session_state.import_jobs[-5:]:  # Show last 5 jobs
            try:
                response = httpx.get(f"{API_BASE}/imaging/import/{job_id}/status")
                if response.status_code == 200:
                    status = response.json()
                    
                    col1, col2, col3 = st.columns([2, 1, 1])
                    
                    with col1:
                        st.text(f"Job {job_id[:8]}...")
                    
                    with col2:
                        progress = status.get("progress_percent", 0)
                        st.progress(progress / 100)
                        st.text(f"{progress:.1f}%")
                    
                    with col3:
                        status_text = status.get("status", "unknown")
                        if status_text == "completed":
                            st.success("✅ Complete")
                        elif status_text == "failed":
                            st.error("❌ Failed")
                        elif status_text == "importing":
                            st.info("🔄 Importing")
                        else:
                            st.warning(f"📋 {status_text}")
                            
            except Exception as e:
                st.error(f"Error checking job {job_id}: {str(e)}")
    else:
        st.info("No recent imports")


def get_available_plates() -> List[str]:
    """Get list of available plates."""
    try:
        # Mock data for now - replace with actual API call
        return ["PLATE_001", "PLATE_002", "PLATE_003"]
    except Exception:
        return []


def get_available_channels(plate_id: Optional[str]) -> List[str]:
    """Get available channels for a plate."""
    try:
        # Mock data for now - replace with actual API call
        return ["DAPI", "GFP", "RFP", "Cy5", "Brightfield"]
    except Exception:
        return []


def get_z_range(plate_id: Optional[str]) -> tuple[int, int]:
    """Get Z-slice range for a plate."""
    try:
        # Mock data for now - replace with actual API call
        return (0, 10)
    except Exception:
        return (0, 0)


def get_t_range(plate_id: Optional[str]) -> tuple[int, int]:
    """Get timepoint range for a plate."""
    try:
        # Mock data for now - replace with actual API call
        return (0, 5)
    except Exception:
        return (0, 0)


def browse_images(
    plate_id: Optional[str],
    well_position: Optional[str],
    channels: List[str],
    z_slice: int,
    timepoint: int,
    qc_filter: str
) -> List[Dict[str, Any]]:
    """Browse images with filters."""
    try:
        params = {
            "z_slice": z_slice,
            "timepoint": timepoint,
            "limit": 50
        }
        
        if plate_id:
            params["plate_id"] = plate_id
        
        if well_position:
            params["well_position"] = well_position
        
        if len(channels) == 1:
            params["channel"] = channels[0]
        
        if qc_filter == "Passed QC only":
            params["passed_qc_only"] = True
        
        response = httpx.get(f"{API_BASE}/imaging/browse", params=params)
        
        if response.status_code == 200:
            data = response.json()
            return data.get("images", [])
        else:
            return []
            
    except Exception as e:
        st.error(f"Error browsing images: {str(e)}")
        return []


def display_image_grid(images: List[Dict[str, Any]]) -> None:
    """Display images in a grid layout."""
    if not images:
        return
    
    # Group by well and channel for organized display
    cols_per_row = 4
    
    for i in range(0, len(images), cols_per_row):
        cols = st.columns(cols_per_row)
        
        for j, col in enumerate(cols):
            if i + j < len(images):
                image = images[i + j]
                
                with col:
                    # Image placeholder (would show actual thumbnail in production)
                    st.image("https://via.placeholder.com/200x200?text=Image", 
                            caption=f"Well {image.get('well_position', 'N/A')} - {image.get('channel', 'N/A')}")
                    
                    # Image metadata
                    st.caption(f"""
                    **Size:** {image.get('width', 0)}x{image.get('height', 0)}
                    **Z:** {image.get('z_slice', 0)} **T:** {image.get('timepoint', 0)}
                    **QC:** {"✅" if image.get('passed_qc') else "❌" if image.get('passed_qc') is False else "❓"}
                    """)
                    
                    if st.button(f"View Details", key=f"view_{image['id']}"):
                        st.session_state.selected_image = image['id']


def get_plate_qc_data(plate_id: str) -> Optional[Dict[str, Any]]:
    """Get QC data for a plate."""
    try:
        response = httpx.get(f"{API_BASE}/imaging/qc/plate/{plate_id}")
        
        if response.status_code == 200:
            return response.json()
        else:
            return None
            
    except Exception as e:
        st.error(f"Error getting QC data: {str(e)}")
        return None


def create_focus_heatmap(focus_data: Dict[str, float]) -> Optional[go.Figure]:
    """Create a focus quality heatmap."""
    if not focus_data:
        return None
    
    try:
        # Parse well positions and create grid
        wells = list(focus_data.keys())
        scores = list(focus_data.values())
        
        # Assume 96-well plate format for now
        rows = 8
        cols = 12
        
        # Create grid matrix
        grid = np.full((rows, cols), np.nan)
        
        for well, score in focus_data.items():
            if len(well) >= 3:  # e.g., "A01"
                row = ord(well[0]) - ord('A')
                col = int(well[1:]) - 1
                
                if 0 <= row < rows and 0 <= col < cols:
                    grid[row, col] = score
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=grid,
            colorscale='RdYlGn',
            zmin=0,
            zmax=1,
            colorbar=dict(title="Focus Score"),
            hoverongaps=False
        ))
        
        fig.update_layout(
            title="Focus Quality by Well Position",
            xaxis_title="Column",
            yaxis_title="Row",
            height=400
        )
        
        # Add well labels
        fig.update_xaxes(
            tickmode='array',
            tickvals=list(range(cols)),
            ticktext=[str(i+1) for i in range(cols)]
        )
        
        fig.update_yaxes(
            tickmode='array',
            tickvals=list(range(rows)),
            ticktext=[chr(ord('A') + i) for i in range(rows)]
        )
        
        return fig
        
    except Exception as e:
        st.error(f"Error creating heatmap: {str(e)}")
        return None


def get_instruments() -> List[Dict[str, Any]]:
    """Get list of registered instruments."""
    try:
        response = httpx.get(f"{API_BASE}/imaging/instruments")
        
        if response.status_code == 200:
            return response.json()
        else:
            return []
            
    except Exception as e:
        st.error(f"Error getting instruments: {str(e)}")
        return []


def register_instrument(name: str, manufacturer: str, model: str, 
                       serial_number: str, facility_location: str) -> None:
    """Register a new instrument."""
    try:
        data = {
            "name": name,
            "manufacturer": manufacturer,
            "model": model
        }
        
        if serial_number:
            data["serial_number"] = serial_number
        
        if facility_location:
            data["facility_location"] = facility_location
        
        response = httpx.post(f"{API_BASE}/imaging/instruments", json=data)
        
        if response.status_code == 201:
            st.success(f"Instrument '{name}' registered successfully!")
            st.rerun()  # Refresh the page
        else:
            st.error(f"Failed to register instrument: {response.text}")
            
    except Exception as e:
        st.error(f"Error registering instrument: {str(e)}")


# Main execution
if __name__ == "__main__":
    render_imaging_browser_page()
