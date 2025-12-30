"""
E2E tests for biophysical assays dashboard page.

This test suite covers the Streamlit dashboard functionality for SPR, MST, and DSC analysis.
"""

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    """Navigate to the Biophysical Assays page."""
    page.goto(f"{base_url}/?page=Biophysical%20Assays")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


@pytest.mark.e2e
def test_page_loads(page: Page, streamlit_server: str):
    """Test that biophysical assays page loads successfully."""
    _goto(page, streamlit_server)
    
    # Check main header is visible
    expect(page.get_by_text("Biophysical Assay Analysis")).to_be_visible(timeout=20000)
    expect(page.get_by_text("Upload SPR, MST, and DSC files")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_tabs_visible(page: Page, streamlit_server: str):
    """Test all 3 tabs are visible and clickable."""
    _goto(page, streamlit_server)
    
    # Check all tabs are present
    expect(page.get_by_text("📈 SPR Analysis")).to_be_visible(timeout=10000)
    expect(page.get_by_text("🌡️ MST Analysis")).to_be_visible(timeout=10000)
    expect(page.get_by_text("🔥 DSC Analysis")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_spr_upload_section(page: Page, streamlit_server: str):
    """Test SPR file upload section exists and functions."""
    _goto(page, streamlit_server)
    
    # Click SPR tab
    page.get_by_text("📈 SPR Analysis").click()
    
    # Check upload components
    expect(page.get_by_text("Surface Plasmon Resonance")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Upload Biacore/Reichert File")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Compound ID (optional)")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Target Name (optional)")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_mst_upload_section(page: Page, streamlit_server: str):
    """Test MST file upload section exists and functions."""
    _goto(page, streamlit_server)
    
    # Click MST tab
    page.get_by_text("🌡️ MST Analysis").click()
    
    # Check upload components
    expect(page.get_by_text("MicroScale Thermophoresis")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Upload NanoTemper File")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_dsc_upload_section(page: Page, streamlit_server: str):
    """Test DSC file upload section exists and functions."""
    _goto(page, streamlit_server)
    
    # Click DSC tab
    page.get_by_text("🔥 DSC Analysis").click()
    
    # Check upload components
    expect(page.get_by_text("Differential Scanning Calorimetry")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Upload MicroCal/TA File")).to_be_visible(timeout=10000)
    expect(page.get_by_text("Protein Name (optional)")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_spr_analysis_section(page: Page, streamlit_server: str):
    """Test SPR analysis section components."""
    _goto(page, streamlit_server)
    
    # Click SPR tab
    page.get_by_text("📈 SPR Analysis").click()
    
    # Check analysis components
    expect(page.get_by_text("Sensorgram Visualization")).to_be_visible(timeout=10000)
    
    # Check for "no experiments" message initially
    expect(page.get_by_text("No SPR experiments found")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_mst_analysis_section(page: Page, streamlit_server: str):
    """Test MST analysis section components."""
    _goto(page, streamlit_server)
    
    # Click MST tab
    page.get_by_text("🌡️ MST Analysis").click()
    
    # Check analysis components
    expect(page.get_by_text("Dose-Response Curve")).to_be_visible(timeout=10000)
    
    # Check for "no experiments" message initially
    expect(page.get_by_text("No MST experiments found")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_dsc_analysis_section(page: Page, streamlit_server: str):
    """Test DSC analysis section components."""
    _goto(page, streamlit_server)
    
    # Click DSC tab
    page.get_by_text("🔥 DSC Analysis").click()
    
    # Check analysis components
    expect(page.get_by_text("Thermogram")).to_be_visible(timeout=10000)
    
    # Check for "no experiments" message initially
    expect(page.get_by_text("No DSC experiments found")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_file_upload_buttons_present(page: Page, streamlit_server: str):
    """Test that file upload buttons are present in all tabs."""
    _goto(page, streamlit_server)
    
    # SPR tab
    page.get_by_text("📈 SPR Analysis").click()
    expect(page.locator("input[type='file']")).to_be_visible(timeout=10000)
    
    # MST tab
    page.get_by_text("🌡️ MST Analysis").click()
    expect(page.locator("input[type='file']")).to_be_visible(timeout=10000)
    
    # DSC tab
    page.get_by_text("🔥 DSC Analysis").click()
    expect(page.locator("input[type='file']")).to_be_visible(timeout=10000)


@pytest.mark.e2e
def test_tab_switching_functionality(page: Page, streamlit_server: str):
    """Test that tab switching works properly."""
    _goto(page, streamlit_server)
    
    # Start with SPR tab
    page.get_by_text("📈 SPR Analysis").click()
    expect(page.get_by_text("Surface Plasmon Resonance")).to_be_visible(timeout=10000)
    
    # Switch to MST tab
    page.get_by_text("🌡️ MST Analysis").click()
    expect(page.get_by_text("MicroScale Thermophoresis")).to_be_visible(timeout=10000)
    
    # Switch to DSC tab
    page.get_by_text("🔥 DSC Analysis").click()
    expect(page.get_by_text("Differential Scanning Calorimetry")).to_be_visible(timeout=10000)
    
    # Switch back to SPR tab
    page.get_by_text("📈 SPR Analysis").click()
    expect(page.get_by_text("Surface Plasmon Resonance")).to_be_visible(timeout=10000)


# Note: Additional tests for file upload functionality, API integration, 
# and data visualization would require mock API responses or test data,
# which are deferred to integration testing with actual backend services.
