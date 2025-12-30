"""Image Analysis Dashboard - High-Content Imaging Workflow."""

from __future__ import annotations

import io
import json
import pandas as pd
import numpy as np
import requests
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List, Optional, Any
from uuid import UUID

from scripts.dashboard.auth import require_auth
from scripts.dashboard.db_session import db_session


def render_image_analysis_page() -> None:
    """Render the Image Analysis page."""
    require_auth()
    
    st.header("🔬 Image Analysis")
    st.caption("High-content imaging workflow for cell morphology and phenotype analysis.")
    
    # Initialize session state
    if "imaging_results" not in st.session_state:
        st.session_state.imaging_results = []
    if "uploaded_images" not in st.session_state:
        st.session_state.uploaded_images = []
    if "selected_plate" not in st.session_state:
        st.session_state.selected_plate = None
    
    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📤 Upload", "🎯 Segment", "📊 Features", "🧪 Plate View"])
    
    with tab1:
        render_upload_tab()
    
    with tab2:
        render_segment_tab()
    
    with tab3:
        render_features_tab()
    
    with tab4:
        render_plate_view_tab()


def render_upload_tab() -> None:
    """Render the image upload tab."""
    st.subheader("Image Upload & Plate Setup")
    st.caption("Upload microscopy images and associate with HTS plates and wells.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Upload Settings")
        
        # Plate information
        plate_barcode = st.text_input(
            "Plate Barcode",
            placeholder="PLATE001",
            help="Unique identifier for the HTS plate"
        )
        
        # Campaign selection
        with db_session() as db:
            try:
                campaigns = _get_hts_campaigns(db)
                campaign_options = {f"{c['campaign_name']} ({c['campaign_id']})": c['id'] for c in campaigns}
                
                if campaign_options:
                    selected_campaign = st.selectbox(
                        "HTS Campaign",
                        options=list(campaign_options.keys()),
                        help="Select the HTS campaign for this plate"
                    )
                    campaign_id = campaign_options.get(selected_campaign)
                else:
                    st.warning("No HTS campaigns found. Create a campaign first.")
                    campaign_id = None
            except Exception as e:
                st.error(f"Failed to load campaigns: {str(e)}")
                campaign_id = None
        
        # Channel selection
        channel = st.selectbox(
            "Channel",
            options=["DAPI", "GFP", "RFP", "Brightfield", "FITC", "TRITC"],
            help="Fluorescence channel or imaging modality"
        )
        
        # Well position
        well_position = st.text_input(
            "Well Position",
            placeholder="A01",
            help="Well position (e.g., A01, B12)"
        )
        
        # Image metadata
        st.markdown("### Image Metadata")
        pixel_size_um = st.number_input(
            "Pixel Size (μm)",
            min_value=0.01,
            max_value=10.0,
            value=0.325,
            step=0.001,
            format="%.3f"
        )
        
        z_slice = st.number_input("Z-slice", min_value=0, value=0)
        timepoint = st.number_input("Timepoint", min_value=0, value=0)
    
    with col2:
        st.markdown("### Image Files")
        
        # File uploader
        uploaded_files = st.file_uploader(
            "Upload Images",
            type=["tiff", "tif", "png", "jpg", "jpeg"],
            accept_multiple_files=True,
            help="Select microscopy image files"
        )
        
        if uploaded_files:
            st.success(f"Uploaded {len(uploaded_files)} files")
            
            # Display file information
            for i, file in enumerate(uploaded_files):
                with st.expander(f"📁 {file.name}"):
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("File Size", f"{file.size / 1024:.1f} KB")
                    with col_b:
                        st.metric("Type", file.type)
            
            # Upload button
            if st.button("🚀 Upload to Platform", type="primary", disabled=not (plate_barcode and campaign_id and well_position)):
                with st.spinner("Uploading images..."):
                    try:
                        results = _upload_images_to_platform(
                            uploaded_files, plate_barcode, campaign_id, 
                            channel, well_position, pixel_size_um, z_slice, timepoint
                        )
                        
                        if results.get("success"):
                            st.success(f"✅ Successfully uploaded {len(uploaded_files)} images")
                            st.session_state.uploaded_images.extend(results.get("image_ids", []))
                            
                            # Show upload summary
                            if results.get("summary"):
                                st.json(results["summary"])
                        else:
                            st.error(f"Upload failed: {results.get('error', 'Unknown error')}")
                    
                    except Exception as e:
                        st.error(f"Upload error: {str(e)}")
        
        # Recent uploads
        if st.session_state.uploaded_images:
            st.markdown("### Recent Uploads")
            for image_id in st.session_state.uploaded_images[-5:]:
                st.text(f"📷 Image ID: {image_id}")


def render_segment_tab() -> None:
    """Render the cell segmentation tab."""
    st.subheader("Cell Segmentation")
    st.caption("Segment cells using CellPose and extract morphological features.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Segmentation Settings")
        
        # Image selection
        if st.session_state.uploaded_images:
            selected_images = st.multiselect(
                "Select Images",
                options=st.session_state.uploaded_images,
                default=st.session_state.uploaded_images[:3] if len(st.session_state.uploaded_images) >= 3 else st.session_state.uploaded_images,
                help="Choose images to segment"
            )
        else:
            st.info("Upload images first to enable segmentation")
            selected_images = []
        
        # Model settings
        model_type = st.selectbox(
            "CellPose Model",
            options=["cyto", "cyto2", "nuclei"],
            index=0,
            help="Pre-trained CellPose model type"
        )
        
        diameter = st.slider(
            "Expected Cell Diameter (pixels)",
            min_value=10,
            max_value=200,
            value=30,
            help="Expected diameter of cells in pixels"
        )
        
        # Channel configuration
        st.markdown("### Channel Configuration")
        cytoplasm_channel = st.selectbox("Cytoplasm Channel", [0, 1, 2], index=0)
        nucleus_channel = st.selectbox("Nucleus Channel", [0, 1, 2], index=0)
        
        # Processing options
        st.markdown("### Processing Options")
        gpu_enabled = st.checkbox("Use GPU", value=True, help="Enable GPU acceleration")
        extract_features = st.checkbox("Extract Features", value=True, help="Calculate morphological features")
        batch_mode = st.checkbox("Batch Processing", value=False, help="Use Celery for background processing")
    
    with col2:
        st.markdown("### Segmentation Preview")
        
        if selected_images:
            # Segmentation controls
            col_a, col_b = st.columns(2)
            
            with col_a:
                if st.button("🎯 Run Segmentation", type="primary", disabled=not selected_images):
                    with st.spinner("Segmenting cells..."):
                        try:
                            if batch_mode and len(selected_images) > 1:
                                # Batch processing
                                result = _run_batch_segmentation(
                                    selected_images, model_type, diameter, 
                                    [cytoplasm_channel, nucleus_channel], extract_features
                                )
                                
                                if result.get("success"):
                                    st.success(f"✅ Batch job queued: {result['task_id']}")
                                    st.info("Check the results tab for progress updates")
                                else:
                                    st.error(f"Batch segmentation failed: {result.get('error')}")
                            
                            else:
                                # Single image processing
                                results = []
                                progress_bar = st.progress(0)
                                
                                for i, image_id in enumerate(selected_images):
                                    result = _run_single_segmentation(
                                        image_id, model_type, diameter,
                                        [cytoplasm_channel, nucleus_channel], extract_features
                                    )
                                    results.append(result)
                                    progress_bar.progress((i + 1) / len(selected_images))
                                
                                # Show results
                                successful = [r for r in results if r.get("success")]
                                st.success(f"✅ Segmented {len(successful)}/{len(selected_images)} images")
                                
                                # Add to session results
                                st.session_state.imaging_results.extend(successful)
                        
                        except Exception as e:
                            st.error(f"Segmentation error: {str(e)}")
            
            with col_b:
                if st.button("🔍 Preview Settings"):
                    st.info("Preview functionality requires image display integration")
        
        # Recent results
        if st.session_state.imaging_results:
            st.markdown("### Recent Segmentations")
            
            results_df = pd.DataFrame([
                {
                    "Image ID": r.get("image_id", "")[:8] + "...",
                    "Cells": r.get("cell_count", 0),
                    "Model": r.get("model_name", ""),
                    "Features": "✅" if r.get("features_extracted") else "❌"
                }
                for r in st.session_state.imaging_results[-10:]
            ])
            
            st.dataframe(results_df, use_container_width=True)


def render_features_tab() -> None:
    """Render the feature analysis tab."""
    st.subheader("Feature Analysis")
    st.caption("Explore extracted morphological features and generate insights.")
    
    # Feature data loading
    if not st.session_state.imaging_results:
        st.info("Run cell segmentation first to analyze features")
        return
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Analysis Settings")
        
        # Segmentation selection
        segmentation_options = {
            f"Seg {i+1} ({r.get('cell_count', 0)} cells)": r.get('segmentation_id')
            for i, r in enumerate(st.session_state.imaging_results)
            if r.get('segmentation_id')
        }
        
        if segmentation_options:
            selected_seg = st.selectbox(
                "Select Segmentation",
                options=list(segmentation_options.keys())
            )
            segmentation_id = segmentation_options[selected_seg]
        else:
            st.warning("No segmentations with extracted features found")
            return
        
        # Feature selection
        available_features = [
            "area", "perimeter", "circularity", 
            "eccentricity", "solidity", "centroid_x", "centroid_y"
        ]
        
        selected_features = st.multiselect(
            "Features to Analyze",
            options=available_features,
            default=["area", "perimeter", "circularity"],
            help="Choose morphological features to visualize"
        )
        
        # Visualization options
        st.markdown("### Visualization Options")
        plot_type = st.radio(
            "Plot Type",
            options=["Histogram", "Scatter Plot", "Box Plot"],
            horizontal=True
        )
        
        # Export options
        st.markdown("### Export")
        if st.button("📊 Export to CSV"):
            try:
                csv_data = _export_features_to_csv(segmentation_id)
                if csv_data:
                    st.download_button(
                        label="Download CSV",
                        data=csv_data,
                        file_name=f"features_{segmentation_id[:8]}.csv",
                        mime="text/csv"
                    )
                else:
                    st.error("No feature data available for export")
            except Exception as e:
                st.error(f"Export failed: {str(e)}")
    
    with col2:
        st.markdown("### Feature Visualizations")
        
        if selected_features and segmentation_id:
            try:
                # Load feature data
                feature_data = _load_feature_data(segmentation_id)
                
                if feature_data and len(feature_data) > 0:
                    df = pd.DataFrame(feature_data)
                    
                    # Generate plots based on selection
                    if plot_type == "Histogram":
                        _render_feature_histograms(df, selected_features)
                    elif plot_type == "Scatter Plot":
                        _render_feature_scatter(df, selected_features)
                    elif plot_type == "Box Plot":
                        _render_feature_boxplots(df, selected_features)
                    
                    # Feature statistics
                    st.markdown("### Feature Statistics")
                    stats_df = df[selected_features].describe()
                    st.dataframe(stats_df, use_container_width=True)
                
                else:
                    st.warning("No feature data found for selected segmentation")
            
            except Exception as e:
                st.error(f"Failed to load feature data: {str(e)}")


def render_plate_view_tab() -> None:
    """Render the plate-level analysis tab."""
    st.subheader("Plate Analysis & Quality Control")
    st.caption("Analyze plate-level patterns, quality metrics, and control performance.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Plate Selection")
        
        # Load available plates
        try:
            with db_session() as db:
                plates = _get_available_plates(db)
                
                if plates:
                    plate_options = {f"{p['barcode']} ({p['well_count']} wells)": p['id'] for p in plates}
                    
                    selected_plate_key = st.selectbox(
                        "Select Plate",
                        options=list(plate_options.keys())
                    )
                    selected_plate_id = plate_options[selected_plate_key]
                    st.session_state.selected_plate = selected_plate_id
                else:
                    st.warning("No plates with imaging data found")
                    return
        except Exception as e:
            st.error(f"Failed to load plates: {str(e)}")
            return
        
        # Feature selection for heatmap
        st.markdown("### Heatmap Settings")
        heatmap_feature = st.selectbox(
            "Feature for Heatmap",
            options=["area_mean", "cell_count", "circularity_mean", "eccentricity_mean"],
            help="Feature to visualize across plate wells"
        )
        
        # Quality control settings
        st.markdown("### Quality Control")
        
        # Control well positions
        positive_controls = st.text_input(
            "Positive Controls",
            placeholder="H01,H02,H23,H24",
            help="Comma-separated list of positive control wells"
        )
        
        negative_controls = st.text_input(
            "Negative Controls",
            placeholder="A01,A02,A23,A24",
            help="Comma-separated list of negative control wells"
        )
        
        # Calculate QC metrics
        if st.button("📊 Calculate QC Metrics"):
            if positive_controls and negative_controls:
                try:
                    pos_wells = [w.strip() for w in positive_controls.split(",")]
                    neg_wells = [w.strip() for w in negative_controls.split(",")]
                    
                    qc_results = _calculate_plate_qc(
                        selected_plate_id, heatmap_feature, pos_wells, neg_wells
                    )
                    
                    if qc_results:
                        st.session_state.qc_results = qc_results
                        st.success("✅ QC metrics calculated")
                    else:
                        st.error("Failed to calculate QC metrics")
                
                except Exception as e:
                    st.error(f"QC calculation failed: {str(e)}")
            else:
                st.warning("Please specify both positive and negative control wells")
    
    with col2:
        st.markdown("### Plate Heatmap")
        
        if st.session_state.selected_plate:
            try:
                # Generate plate heatmap
                heatmap_data = _get_plate_heatmap_data(st.session_state.selected_plate, heatmap_feature)
                
                if heatmap_data:
                    _render_plate_heatmap(heatmap_data, heatmap_feature)
                else:
                    st.info("No data available for plate heatmap")
            
            except Exception as e:
                st.error(f"Failed to generate heatmap: {str(e)}")
        
        # QC Results
        if hasattr(st.session_state, 'qc_results') and st.session_state.qc_results:
            st.markdown("### Quality Control Results")
            
            qc = st.session_state.qc_results
            
            # Z' factor display
            if 'zprime' in qc and qc['zprime']:
                zprime_val = qc['zprime'].get('zprime_factor', 0)
                quality = qc['zprime'].get('assay_quality', 'Unknown')
                
                col_a, col_b = st.columns(2)
                with col_a:
                    st.metric("Z' Factor", f"{zprime_val:.3f}")
                with col_b:
                    st.metric("Assay Quality", quality)
                
                # Color-code quality
                if zprime_val > 0.5:
                    st.success("Excellent assay quality (Z' > 0.5)")
                elif zprime_val > 0:
                    st.warning("Marginal assay quality (0 < Z' < 0.5)")
                else:
                    st.error("Poor assay quality (Z' < 0)")
            
            # Uniformity metrics
            if 'uniformity' in qc and qc['uniformity']:
                uniformity = qc['uniformity']
                
                st.markdown("#### Plate Uniformity")
                col_a, col_b, col_c = st.columns(3)
                
                with col_a:
                    st.metric("CV (%)", f"{uniformity.get('cv_percent', 0):.1f}")
                with col_b:
                    st.metric("Mean", f"{uniformity.get('mean', 0):.1f}")
                with col_c:
                    st.metric("Std Dev", f"{uniformity.get('std', 0):.1f}")
            
            # Outliers
            if 'outliers_2sigma' in qc:
                outliers = qc['outliers_2sigma']
                if outliers:
                    st.warning(f"Outlier wells (2σ): {', '.join(outliers)}")
                else:
                    st.success("No outlier wells detected")


# Helper functions

def _get_hts_campaigns(db) -> List[Dict]:
    """Get available HTS campaigns."""
    try:
        from amprenta_rag.models.chemistry import HTSCampaign
        campaigns = db.query(HTSCampaign).all()
        return [
            {
                "id": str(c.id),
                "campaign_id": c.campaign_id,
                "campaign_name": c.campaign_name,
                "description": c.description
            }
            for c in campaigns
        ]
    except Exception:
        return []


def _upload_images_to_platform(files, plate_barcode, campaign_id, channel, well_position, pixel_size_um, z_slice, timepoint) -> Dict:
    """Upload images to the platform via API."""
    try:
        # Mock implementation - would call actual API
        return {
            "success": True,
            "image_ids": [f"img_{i}_{file.name[:8]}" for i, file in enumerate(files)],
            "summary": {
                "plate_barcode": plate_barcode,
                "files_uploaded": len(files),
                "channel": channel,
                "well_position": well_position
            }
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


def _run_batch_segmentation(image_ids, model_type, diameter, channels, extract_features) -> Dict:
    """Run batch segmentation via Celery."""
    try:
        # Mock implementation - would call /api/v1/imaging/segment-batch
        return {
            "success": True,
            "task_id": f"task_{len(image_ids)}_{model_type}",
            "image_count": len(image_ids)
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


def _run_single_segmentation(image_id, model_type, diameter, channels, extract_features) -> Dict:
    """Run single image segmentation."""
    try:
        # Mock implementation - would call /api/v1/imaging/segment
        import random
        return {
            "success": True,
            "image_id": image_id,
            "segmentation_id": f"seg_{image_id[:8]}",
            "cell_count": random.randint(50, 200),
            "model_name": model_type,
            "features_extracted": extract_features
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


def _load_feature_data(segmentation_id) -> List[Dict]:
    """Load feature data for a segmentation."""
    # Mock implementation - would call /api/v1/imaging/features/{segmentation_id}
    import random
    np.random.seed(42)
    
    n_cells = random.randint(50, 150)
    return [
        {
            "cell_id": i + 1,
            "area": np.random.normal(100, 20),
            "perimeter": np.random.normal(40, 8),
            "circularity": np.random.beta(2, 2),
            "eccentricity": np.random.beta(2, 5),
            "solidity": np.random.beta(5, 2),
            "centroid_x": np.random.uniform(0, 1024),
            "centroid_y": np.random.uniform(0, 1024)
        }
        for i in range(n_cells)
    ]


def _export_features_to_csv(segmentation_id) -> Optional[str]:
    """Export feature data to CSV format."""
    try:
        feature_data = _load_feature_data(segmentation_id)
        if feature_data:
            df = pd.DataFrame(feature_data)
            return df.to_csv(index=False)
        return None
    except Exception:
        return None


def _get_available_plates(db) -> List[Dict]:
    """Get plates with imaging data."""
    try:
        # Mock implementation - would query actual database
        return [
            {
                "id": "plate_001",
                "barcode": "PLATE001",
                "well_count": 384,
                "image_count": 1152
            },
            {
                "id": "plate_002", 
                "barcode": "PLATE002",
                "well_count": 96,
                "image_count": 288
            }
        ]
    except Exception:
        return []


def _get_plate_heatmap_data(plate_id, feature) -> Optional[Dict]:
    """Get heatmap data for plate visualization."""
    # Mock implementation - would call /api/v1/imaging/plates/{plate_id}/heatmap
    try:
        np.random.seed(42)
        
        # Generate mock 96-well plate data
        wells = []
        for row in "ABCDEFGH":
            for col in range(1, 13):
                wells.append(f"{row}{col:02d}")
        
        return {
            well: np.random.normal(100, 20) if np.random.random() > 0.1 else np.nan
            for well in wells
        }
    except Exception:
        return None


def _calculate_plate_qc(plate_id, feature, pos_controls, neg_controls) -> Optional[Dict]:
    """Calculate plate quality control metrics."""
    # Mock implementation - would call QC API endpoints
    try:
        import random
        
        # Mock Z' factor calculation
        zprime = random.uniform(-0.5, 0.8)
        quality = "Excellent" if zprime > 0.5 else "Marginal" if zprime > 0 else "Poor"
        
        return {
            "zprime": {
                "zprime_factor": zprime,
                "assay_quality": quality,
                "positive_controls": {"mean": 150, "std": 10},
                "negative_controls": {"mean": 50, "std": 8}
            },
            "uniformity": {
                "cv_percent": random.uniform(5, 25),
                "mean": random.uniform(80, 120),
                "std": random.uniform(10, 30)
            },
            "outliers_2sigma": random.sample(["A01", "B05", "G12"], random.randint(0, 2))
        }
    except Exception:
        return None


def _render_feature_histograms(df, features):
    """Render feature histograms."""
    for feature in features[:4]:  # Limit to 4 features
        if feature in df.columns:
            fig = px.histogram(
                df, x=feature, 
                title=f"Distribution of {feature.title()}",
                nbins=30
            )
            st.plotly_chart(fig, use_container_width=True)


def _render_feature_scatter(df, features):
    """Render feature scatter plots."""
    if len(features) >= 2:
        fig = px.scatter(
            df, x=features[0], y=features[1],
            title=f"{features[0].title()} vs {features[1].title()}"
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Select at least 2 features for scatter plot")


def _render_feature_boxplots(df, features):
    """Render feature box plots."""
    if features:
        fig = go.Figure()
        for feature in features[:5]:  # Limit to 5 features
            if feature in df.columns:
                fig.add_trace(go.Box(y=df[feature], name=feature.title()))
        
        fig.update_layout(title="Feature Distributions")
        st.plotly_chart(fig, use_container_width=True)


def _render_plate_heatmap(heatmap_data, feature):
    """Render plate heatmap visualization."""
    try:
        # Convert well positions to coordinates
        wells = []
        values = []
        
        for well, value in heatmap_data.items():
            if pd.notna(value):
                row = ord(well[0]) - ord('A')
                col = int(well[1:]) - 1
                wells.append((row, col, well))
                values.append(value)
        
        if not wells:
            st.warning("No data available for heatmap")
            return
        
        # Create heatmap matrix
        max_row = max(w[0] for w in wells) + 1
        max_col = max(w[1] for w in wells) + 1
        
        heatmap_matrix = np.full((max_row, max_col), np.nan)
        well_labels = [[""] * max_col for _ in range(max_row)]
        
        for (row, col, well), value in zip(wells, values):
            heatmap_matrix[row, col] = value
            well_labels[row][col] = well
        
        # Create plotly heatmap
        fig = px.imshow(
            heatmap_matrix,
            title=f"Plate Heatmap - {feature.title()}",
            color_continuous_scale="viridis",
            aspect="auto"
        )
        
        # Add well labels
        fig.update_traces(
            text=well_labels,
            texttemplate="%{text}",
            textfont={"size": 8}
        )
        
        fig.update_layout(
            xaxis_title="Column",
            yaxis_title="Row"
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
    except Exception as e:
        st.error(f"Failed to render heatmap: {str(e)}")


if __name__ == "__main__":
    render_image_analysis_page()
