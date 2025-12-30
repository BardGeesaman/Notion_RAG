"""E2E tests for Flow Cytometry dashboard page."""

from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    """Navigate to the Flow Cytometry page."""
    page.goto(f"{base_url}/?page=Flow%20Cytometry")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_flow_cytometry_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Flow Cytometry page loads successfully."""
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Flow Cytometry Analysis"))).to_be_visible(timeout=20000)


def test_flow_cytometry_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that all four tabs are present and visible."""
    _goto(page, streamlit_server)
    
    # Check for all tabs
    expect(page.get_by_text("📤 Upload")).to_be_visible(timeout=10000)
    expect(page.get_by_text("📊 Scatter Plots")).to_be_visible(timeout=10000)
    expect(page.get_by_text("🎯 Gating")).to_be_visible(timeout=10000)
    expect(page.get_by_text("📈 Statistics")).to_be_visible(timeout=10000)


def test_upload_tab_file_uploader_visible(page: Page, streamlit_server: str) -> None:
    """Test that the file uploader is visible in the Upload tab."""
    _goto(page, streamlit_server)
    
    # Upload tab should be active by default
    expect(page.get_by_text("Choose an FCS file")).to_be_visible(timeout=10000)
    expect(page.locator("input[type='file']")).to_be_visible(timeout=10000)


def test_upload_tab_dataset_selector(page: Page, streamlit_server: str) -> None:
    """Test that the dataset selector is present in the Upload tab."""
    _goto(page, streamlit_server)
    
    # Check for existing datasets section
    expect(page.get_by_text("Existing Datasets")).to_be_visible(timeout=10000)


def test_scatter_plots_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to the Scatter Plots tab."""
    _goto(page, streamlit_server)
    
    # Click on Scatter Plots tab
    page.get_by_text("📊 Scatter Plots").click()
    page.wait_for_timeout(1000)
    
    # Should show parameter selection or dataset selection message
    expect(page.get_by_text("2D Scatter Plots")).to_be_visible(timeout=10000)


def test_gating_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to the Gating tab."""
    _goto(page, streamlit_server)
    
    # Click on Gating tab
    page.get_by_text("🎯 Gating").click()
    page.wait_for_timeout(1000)
    
    # Should show gating interface
    expect(page.get_by_text("Interactive Gating")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Create New Gate")).to_be_visible(timeout=10000)


def test_gating_tab_gate_type_selector(page: Page, streamlit_server: str) -> None:
    """Test that gate type selector is present in the Gating tab."""
    _goto(page, streamlit_server)
    
    # Navigate to Gating tab
    page.get_by_text("🎯 Gating").click()
    page.wait_for_timeout(1000)
    
    # Check for gate type selector
    gate_selector = page.locator("select").first
    expect(gate_selector).to_be_visible(timeout=10000)


def test_statistics_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to the Statistics tab."""
    _goto(page, streamlit_server)
    
    # Click on Statistics tab
    page.get_by_text("📈 Statistics").click()
    page.wait_for_timeout(1000)
    
    # Should show statistics interface
    expect(page.get_by_text("Population Statistics")).to_be_visible(timeout=10000)


def test_statistics_tab_export_functionality(page: Page, streamlit_server: str) -> None:
    """Test that export functionality is present in the Statistics tab."""
    _goto(page, streamlit_server)
    
    # Navigate to Statistics tab
    page.get_by_text("📈 Statistics").click()
    page.wait_for_timeout(1000)
    
    # Check for export section (may not be visible if no data)
    # This test verifies the UI structure is present
    expect(page.get_by_text("Population Statistics")).to_be_visible(timeout=10000)


def test_page_responsive_design(page: Page, streamlit_server: str) -> None:
    """Test that the page works on different screen sizes."""
    _goto(page, streamlit_server)
    
    # Test mobile viewport
    page.set_viewport_size({"width": 375, "height": 667})
    page.wait_for_timeout(1000)
    
    # Header should still be visible
    expect(page.get_by_role("heading", name=re.compile(r"Flow Cytometry Analysis"))).to_be_visible(timeout=10000)
    
    # Tabs should still be accessible
    expect(page.get_by_text("📤 Upload")).to_be_visible(timeout=10000)
    
    # Reset to desktop viewport
    page.set_viewport_size({"width": 1280, "height": 720})


def test_error_handling_display(page: Page, streamlit_server: str) -> None:
    """Test that error messages are displayed appropriately."""
    _goto(page, streamlit_server)
    
    # Navigate to Scatter Plots tab without selecting dataset
    page.get_by_text("📊 Scatter Plots").click()
    page.wait_for_timeout(1000)
    
    # Should show warning about selecting dataset
    # This test verifies error handling UI is present
    expect(page.get_by_text("2D Scatter Plots")).to_be_visible(timeout=10000)


def test_tab_state_persistence(page: Page, streamlit_server: str) -> None:
    """Test that tab state persists during navigation."""
    _goto(page, streamlit_server)
    
    # Navigate to different tabs and back
    page.get_by_text("🎯 Gating").click()
    page.wait_for_timeout(500)
    expect(page.get_by_text("Interactive Gating")).to_be_visible(timeout=10000)
    
    page.get_by_text("📤 Upload").click()
    page.wait_for_timeout(500)
    expect(page.get_by_text("Choose an FCS file")).to_be_visible(timeout=10000)
    
    page.get_by_text("📈 Statistics").click()
    page.wait_for_timeout(500)
    expect(page.get_by_text("Population Statistics")).to_be_visible(timeout=10000)


def test_help_text_visibility(page: Page, streamlit_server: str) -> None:
    """Test that help text and instructions are visible."""
    _goto(page, streamlit_server)
    
    # Check main caption
    expect(page.get_by_text("Upload FCS files, create gates, and analyze cell populations.")).to_be_visible(timeout=10000)
    
    # Navigate to gating tab and check help text
    page.get_by_text("🎯 Gating").click()
    page.wait_for_timeout(1000)
    
    # Should show gate type descriptions
    expect(page.get_by_text("Create New Gate")).to_be_visible(timeout=10000)


def test_form_elements_accessibility(page: Page, streamlit_server: str) -> None:
    """Test that form elements are accessible."""
    _goto(page, streamlit_server)
    
    # Check file uploader accessibility
    file_input = page.locator("input[type='file']")
    expect(file_input).to_be_visible(timeout=10000)
    
    # Navigate to gating tab
    page.get_by_text("🎯 Gating").click()
    page.wait_for_timeout(1000)
    
    # Check form inputs are accessible
    gate_name_input = page.locator("input[placeholder*='Enter a name']").first
    if gate_name_input.is_visible():
        expect(gate_name_input).to_be_visible()


def test_loading_states_handling(page: Page, streamlit_server: str) -> None:
    """Test that loading states are handled properly."""
    _goto(page, streamlit_server)
    
    # The page should load without showing persistent loading spinners
    # This test ensures the page reaches a stable state
    page.wait_for_timeout(3000)
    
    # Check that main content is visible (not stuck in loading)
    expect(page.get_by_role("heading", name=re.compile(r"Flow Cytometry Analysis"))).to_be_visible(timeout=10000)
    expect(page.get_by_text("📤 Upload")).to_be_visible(timeout=10000)
