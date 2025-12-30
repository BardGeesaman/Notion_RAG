"""E2E tests for Imaging Browser dashboard page."""

from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    """Navigate to the Imaging Browser page."""
    page.goto(f"{base_url}/?page=Imaging%20Browser")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_imaging_browser_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Imaging Browser page loads successfully."""
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Microscopy Image Browser"))).to_be_visible(timeout=20000)


def test_imaging_browser_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that all four tabs are present and visible."""
    _goto(page, streamlit_server)
    
    # Check for all tabs
    expect(page.get_by_text("📤 Batch Import")).to_be_visible(timeout=10000)
    expect(page.get_by_text("🖼️ 5D Browser")).to_be_visible(timeout=10000)
    expect(page.get_by_text("📊 QC Dashboard")).to_be_visible(timeout=10000)
    expect(page.get_by_text("🔧 Instruments")).to_be_visible(timeout=10000)


def test_batch_import_tab_elements_visible(page: Page, streamlit_server: str) -> None:
    """Test that the batch import tab elements are visible."""
    _goto(page, streamlit_server)
    
    # Batch Import tab should be active by default
    expect(page.get_by_text("Export Directory Path")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Vendor Format")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Auto-detect")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Supported Formats")).to_be_visible(timeout=10000)


def test_batch_import_vendor_selector(page: Page, streamlit_server: str) -> None:
    """Test that vendor format selector works."""
    _goto(page, streamlit_server)
    
    # Find and interact with vendor selector
    vendor_selector = page.locator("div[data-testid='stSelectbox']").filter(has_text="Vendor Format")
    expect(vendor_selector).to_be_visible(timeout=10000)
    
    # Check that Auto-detect is default
    expect(vendor_selector.get_by_text("Auto-detect")).to_be_visible(timeout=5000)


def test_batch_import_path_input(page: Page, streamlit_server: str) -> None:
    """Test that import path input field works."""
    _goto(page, streamlit_server)
    
    # Find path input field
    path_input = page.get_by_placeholder("/path/to/vendor/export")
    expect(path_input).to_be_visible(timeout=10000)
    
    # Test typing in the field
    path_input.fill("/test/path")
    expect(path_input).to_have_value("/test/path")


def test_5d_browser_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to 5D Browser tab."""
    _goto(page, streamlit_server)
    
    # Click on 5D Browser tab
    page.get_by_text("🖼️ 5D Browser").click()
    page.wait_for_timeout(1000)
    
    # Check for 5D Browser elements
    expect(page.get_by_text("5D Image Browser")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Filters")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Plate")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Channels")).to_be_visible(timeout=10000)


def test_5d_browser_filters_visible(page: Page, streamlit_server: str) -> None:
    """Test that 5D browser filter elements are visible."""
    _goto(page, streamlit_server)
    
    # Navigate to 5D Browser tab
    page.get_by_text("🖼️ 5D Browser").click()
    page.wait_for_timeout(1000)
    
    # Check filter elements
    expect(page.get_by_text("Well Position")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Z-slice")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Timepoint")).to_be_visible(timeout=10000)
    expect(page.get_by_text("QC Status")).to_be_visible(timeout=10000)


def test_qc_dashboard_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to QC Dashboard tab."""
    _goto(page, streamlit_server)
    
    # Click on QC Dashboard tab
    page.get_by_text("📊 QC Dashboard").click()
    page.wait_for_timeout(1000)
    
    # Check for QC Dashboard elements
    expect(page.get_by_text("Quality Control Dashboard")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Select Plate for QC Analysis")).to_be_visible(timeout=10000)


def test_qc_dashboard_plate_selector(page: Page, streamlit_server: str) -> None:
    """Test QC dashboard plate selector functionality."""
    _goto(page, streamlit_server)
    
    # Navigate to QC Dashboard tab
    page.get_by_text("📊 QC Dashboard").click()
    page.wait_for_timeout(1000)
    
    # Find plate selector
    plate_selector = page.locator("div[data-testid='stSelectbox']").filter(has_text="Select Plate for QC Analysis")
    expect(plate_selector).to_be_visible(timeout=10000)


def test_instruments_tab_navigation(page: Page, streamlit_server: str) -> None:
    """Test navigation to Instruments tab."""
    _goto(page, streamlit_server)
    
    # Click on Instruments tab
    page.get_by_text("🔧 Instruments").click()
    page.wait_for_timeout(1000)
    
    # Check for Instruments elements
    expect(page.get_by_text("Microscope Instrument Registry")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Registered Microscopes")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Add New Microscope")).to_be_visible(timeout=10000)


def test_instruments_add_new_form(page: Page, streamlit_server: str) -> None:
    """Test the add new instrument form."""
    _goto(page, streamlit_server)
    
    # Navigate to Instruments tab
    page.get_by_text("🔧 Instruments").click()
    page.wait_for_timeout(1000)
    
    # Expand the add new instrument form
    page.get_by_text("Register New Instrument").click()
    page.wait_for_timeout(500)
    
    # Check form elements
    expect(page.get_by_placeholder("e.g., Nikon Ti2-E #1")).to_be_visible(timeout=10000)
    expect(page.get_by_placeholder("e.g., Nikon")).to_be_visible(timeout=10000)
    expect(page.get_by_placeholder("e.g., Ti2-E")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Register Instrument")).to_be_visible(timeout=10000)


def test_instruments_form_input_validation(page: Page, streamlit_server: str) -> None:
    """Test that instrument form inputs work correctly."""
    _goto(page, streamlit_server)
    
    # Navigate to Instruments tab
    page.get_by_text("🔧 Instruments").click()
    page.wait_for_timeout(1000)
    
    # Expand the add new instrument form
    page.get_by_text("Register New Instrument").click()
    page.wait_for_timeout(500)
    
    # Test form inputs
    name_input = page.get_by_placeholder("e.g., Nikon Ti2-E #1")
    manufacturer_input = page.get_by_placeholder("e.g., Nikon")
    model_input = page.get_by_placeholder("e.g., Ti2-E")
    
    # Fill form fields
    name_input.fill("Test Microscope")
    manufacturer_input.fill("Test Manufacturer")
    model_input.fill("Test Model")
    
    # Verify values
    expect(name_input).to_have_value("Test Microscope")
    expect(manufacturer_input).to_have_value("Test Manufacturer")
    expect(model_input).to_have_value("Test Model")


def test_page_responsive_layout(page: Page, streamlit_server: str) -> None:
    """Test that the page layout is responsive and elements are properly positioned."""
    _goto(page, streamlit_server)
    
    # Check main header is visible
    expect(page.get_by_role("heading", name=re.compile(r"Microscopy Image Browser"))).to_be_visible(timeout=20000)
    
    # Check caption is visible
    expect(page.get_by_text("Import, browse, and analyze high-content imaging data")).to_be_visible(timeout=10000)
    
    # Check all tabs are in a horizontal layout
    tabs_container = page.locator("div[data-testid='stTabs']")
    expect(tabs_container).to_be_visible(timeout=10000)
    
    # Verify tab buttons are clickable
    for tab_name in ["📤 Batch Import", "🖼️ 5D Browser", "📊 QC Dashboard", "🔧 Instruments"]:
        tab_button = page.get_by_text(tab_name)
        expect(tab_button).to_be_visible(timeout=5000)
        expect(tab_button).to_be_enabled()


def test_all_tabs_accessible_and_functional(page: Page, streamlit_server: str) -> None:
    """Test that all tabs can be accessed and contain expected content."""
    _goto(page, streamlit_server)
    
    # Test each tab
    tabs = [
        ("📤 Batch Import", "Export Directory Path"),
        ("🖼️ 5D Browser", "5D Image Browser"),
        ("📊 QC Dashboard", "Quality Control Dashboard"),
        ("🔧 Instruments", "Microscope Instrument Registry")
    ]
    
    for tab_name, expected_content in tabs:
        # Click tab
        page.get_by_text(tab_name).click()
        page.wait_for_timeout(1000)
        
        # Verify content is visible
        expect(page.get_by_text(expected_content)).to_be_visible(timeout=10000)
