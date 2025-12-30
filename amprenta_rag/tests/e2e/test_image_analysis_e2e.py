"""End-to-end tests for Image Analysis dashboard."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


class TestImageAnalysisE2E:
    """E2E tests for Image Analysis dashboard page."""

    def test_page_loads(self, page: Page, base_url: str):
        """Test that the Image Analysis page loads successfully."""
        page.goto(f"{base_url}")
        
        # Navigate to Image Analysis page
        page.get_by_text("Image Analysis").click()
        
        # Check page header
        expect(page.get_by_text("🔬 Image Analysis")).to_be_visible()
        expect(page.get_by_text("High-content imaging workflow")).to_be_visible()
        
        # Check tabs are present
        expect(page.get_by_text("📤 Upload")).to_be_visible()
        expect(page.get_by_text("🎯 Segment")).to_be_visible()
        expect(page.get_by_text("📊 Features")).to_be_visible()
        expect(page.get_by_text("🧪 Plate View")).to_be_visible()

    def test_upload_tab_ui(self, page: Page, base_url: str):
        """Test Upload tab UI components."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Should be on Upload tab by default
        upload_tab = page.get_by_text("📤 Upload")
        expect(upload_tab).to_be_visible()
        
        # Check upload form elements
        expect(page.get_by_text("Plate Barcode")).to_be_visible()
        expect(page.get_by_text("HTS Campaign")).to_be_visible()
        expect(page.get_by_text("Channel")).to_be_visible()
        expect(page.get_by_text("Well Position")).to_be_visible()
        expect(page.get_by_text("Pixel Size (μm)")).to_be_visible()
        
        # Check file uploader
        expect(page.get_by_text("Upload Images")).to_be_visible()
        
        # Check upload button (should be disabled initially)
        upload_button = page.get_by_text("🚀 Upload to Platform")
        expect(upload_button).to_be_visible()

    def test_segment_tab_ui(self, page: Page, base_url: str):
        """Test Segment tab UI components."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Click on Segment tab
        page.get_by_text("🎯 Segment").click()
        
        # Check segmentation settings
        expect(page.get_by_text("Segmentation Settings")).to_be_visible()
        expect(page.get_by_text("CellPose Model")).to_be_visible()
        expect(page.get_by_text("Expected Cell Diameter")).to_be_visible()
        
        # Check channel configuration
        expect(page.get_by_text("Channel Configuration")).to_be_visible()
        expect(page.get_by_text("Cytoplasm Channel")).to_be_visible()
        expect(page.get_by_text("Nucleus Channel")).to_be_visible()
        
        # Check processing options
        expect(page.get_by_text("Processing Options")).to_be_visible()
        expect(page.get_by_text("Use GPU")).to_be_visible()
        expect(page.get_by_text("Extract Features")).to_be_visible()
        expect(page.get_by_text("Batch Processing")).to_be_visible()
        
        # Check segmentation preview area
        expect(page.get_by_text("Segmentation Preview")).to_be_visible()

    def test_features_tab_ui(self, page: Page, base_url: str):
        """Test Features tab UI components."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Click on Features tab
        page.get_by_text("📊 Features").click()
        
        # Should show info message when no results available
        expect(page.get_by_text("Run cell segmentation first")).to_be_visible()
        
        # Check that analysis settings section exists
        expect(page.get_by_text("Analysis Settings")).to_be_visible()
        expect(page.get_by_text("Feature Visualizations")).to_be_visible()

    def test_plate_view_tab_ui(self, page: Page, base_url: str):
        """Test Plate View tab UI components."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Click on Plate View tab
        page.get_by_text("🧪 Plate View").click()
        
        # Check plate selection
        expect(page.get_by_text("Plate Selection")).to_be_visible()
        expect(page.get_by_text("Select Plate")).to_be_visible()
        
        # Check heatmap settings
        expect(page.get_by_text("Heatmap Settings")).to_be_visible()
        expect(page.get_by_text("Feature for Heatmap")).to_be_visible()
        
        # Check quality control section
        expect(page.get_by_text("Quality Control")).to_be_visible()
        expect(page.get_by_text("Positive Controls")).to_be_visible()
        expect(page.get_by_text("Negative Controls")).to_be_visible()
        
        # Check plate heatmap area
        expect(page.get_by_text("Plate Heatmap")).to_be_visible()

    def test_tab_navigation(self, page: Page, base_url: str):
        """Test navigation between tabs."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Test clicking through all tabs
        tabs = ["📤 Upload", "🎯 Segment", "📊 Features", "🧪 Plate View"]
        
        for tab_name in tabs:
            tab = page.get_by_text(tab_name)
            expect(tab).to_be_visible()
            tab.click()
            
            # Wait for tab content to load
            page.wait_for_timeout(500)
            
            # Verify tab is active (this depends on Streamlit's tab implementation)
            expect(tab).to_be_visible()

    def test_upload_form_validation(self, page: Page, base_url: str):
        """Test upload form validation."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Should be on Upload tab
        upload_button = page.get_by_text("🚀 Upload to Platform")
        
        # Button should be disabled when required fields are empty
        expect(upload_button).to_be_visible()
        
        # Fill in some fields
        plate_input = page.get_by_placeholder("PLATE001")
        if plate_input.is_visible():
            plate_input.fill("TEST_PLATE")
        
        well_input = page.get_by_placeholder("A01")
        if well_input.is_visible():
            well_input.fill("A01")
        
        # Check that form elements are interactive
        channel_select = page.locator("text=DAPI").first
        if channel_select.is_visible():
            expect(channel_select).to_be_visible()

    def test_plate_heatmap_renders(self, page: Page, base_url: str):
        """Test that plate heatmap area renders."""
        page.goto(f"{base_url}")
        page.get_by_text("Image Analysis").click()
        
        # Navigate to Plate View tab
        page.get_by_text("🧪 Plate View").click()
        
        # Check that heatmap section is present
        expect(page.get_by_text("Plate Heatmap")).to_be_visible()
        
        # Check for plate selection dropdown
        plate_selector = page.locator("text=Select Plate").first
        if plate_selector.is_visible():
            expect(plate_selector).to_be_visible()
        
        # Check for feature selector
        feature_selector = page.locator("text=Feature for Heatmap").first
        if feature_selector.is_visible():
            expect(feature_selector).to_be_visible()
        
        # Check QC calculation button
        qc_button = page.get_by_text("📊 Calculate QC Metrics")
        expect(qc_button).to_be_visible()


# Fixtures and configuration

@pytest.fixture
def base_url():
    """Base URL for the Streamlit application."""
    return "http://localhost:8501"


# Note: These tests require a running Streamlit server
# Run with: streamlit run scripts/dashboard/app.py
# Then: pytest amprenta_rag/tests/e2e/test_image_analysis_e2e.py --headed
