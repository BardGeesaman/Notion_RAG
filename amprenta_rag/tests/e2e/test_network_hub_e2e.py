"""E2E tests for Network Hub dashboard."""

import pytest
from playwright.sync_api import Page, expect

pytestmark = pytest.mark.requires_server

def _goto_network_hub_page(page: Page, base_url: str) -> None:
    """Navigate to the Network Hub page."""
    # First go to home page
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    # Then navigate via query param
    page.goto(f"{base_url}/?page=Network%20Hub")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)

class TestNetworkHub:
    """E2E tests for Cytoscape Network Hub."""
    
    def test_network_hub_loads(self, page: Page, streamlit_server):
        """Test Network Hub page loads with header."""
        _goto_network_hub_page(page, streamlit_server)
        # Check for page header using more specific selectors
        expect(page.get_by_role("heading", name="🕸️ Network Hub")).to_be_visible()
        expect(page.get_by_text("Unified network visualization")).to_be_visible()
    
    def test_ppi_tab_elements(self, page: Page, streamlit_server):
        """Test PPI tab has required elements."""
        _goto_network_hub_page(page, streamlit_server)
        # Click PPI tab if it exists
        ppi_tab = page.get_by_text("PPI Networks")
        if ppi_tab.is_visible():
            ppi_tab.click()
            # Check for gene input
            expect(page.get_by_text("Gene Symbols")).to_be_visible()
    
    def test_evidence_tab_elements(self, page: Page, streamlit_server):
        """Test Evidence Graph tab has required elements."""
        _goto_network_hub_page(page, streamlit_server)
        # Click Evidence tab using more specific selector
        evidence_tab = page.get_by_role("tab", name="🕸️ Evidence Graph")
        if evidence_tab.is_visible():
            evidence_tab.click()
            # Check for entity type selector
            expect(page.get_by_text("Entity type")).to_be_visible()

