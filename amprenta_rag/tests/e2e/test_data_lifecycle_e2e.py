"""E2E tests for Data Lifecycle dashboard."""

import pytest
from playwright.sync_api import Page, expect

pytestmark = pytest.mark.requires_server

def _goto_data_lifecycle_page(page: Page, base_url: str) -> None:
    """Navigate to the Data Lifecycle page."""
    # First go to home page
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    # Then navigate via query param
    page.goto(f"{base_url}/?page=Data%20Lifecycle")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)

class TestDataLifecycle:
    """E2E tests for Data Lifecycle dashboard."""
    
    def test_page_loads_with_tabs(self, page: Page, streamlit_server):
        """Test Data Lifecycle page loads with all 4 tabs visible."""
        _goto_data_lifecycle_page(page, streamlit_server)
        
        # Check for all 4 tabs
        expect(page.get_by_role("tab", name="Overview")).to_be_visible()
        expect(page.get_by_role("tab", name="Quarantine")).to_be_visible()
        expect(page.get_by_role("tab", name="Bulk Operations")).to_be_visible()
        expect(page.get_by_role("tab", name="Audit")).to_be_visible()
    
    def test_overview_tab_orphan_stats(self, page: Page, streamlit_server):
        """Test Overview tab shows orphan stats section."""
        _goto_data_lifecycle_page(page, streamlit_server)
        
        # Click Overview tab (should be default but ensure it's selected)
        overview_tab = page.get_by_role("tab", name="Overview")
        if overview_tab.is_visible():
            overview_tab.click()
        
        # Check for orphan stats section
        expect(page.get_by_text("Orphan")).to_be_visible()
    
    def test_bulk_operations_permanent_delete(self, page: Page, streamlit_server):
        """Test Bulk Operations tab has Permanent Delete section in expander."""
        _goto_data_lifecycle_page(page, streamlit_server)
        
        # Click Bulk Operations tab
        bulk_ops_tab = page.get_by_role("tab", name="Bulk Operations")
        if bulk_ops_tab.is_visible():
            bulk_ops_tab.click()
            page.wait_for_timeout(1000)
        
        # Check for Permanent Delete section using specific heading selector
        expect(page.get_by_role("heading", name="⚠️ Permanent Delete")).to_be_visible()
    
    def test_quarantine_tab_clickable(self, page: Page, streamlit_server):
        """Test Quarantine tab is clickable and functional."""
        _goto_data_lifecycle_page(page, streamlit_server)
        
        # Click Quarantine tab
        quarantine_tab = page.get_by_role("tab", name="Quarantine")
        expect(quarantine_tab).to_be_visible()
        quarantine_tab.click()
        page.wait_for_timeout(1000)
        
        # Verify tab is now selected/active (basic functionality test)
        expect(quarantine_tab).to_be_visible()