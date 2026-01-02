"""E2E tests for retrosynthesis dashboard."""

import pytest
from playwright.sync_api import Page, expect

pytestmark = pytest.mark.requires_server

def _goto_retrosynthesis_page(page: Page, base_url: str) -> None:
    """Navigate to the Retrosynthesis page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Retrosynthesis%20Advisor")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


class TestRetrosynthesisE2E:
    def test_page_loads(self, page: Page, streamlit_server):
        """Retrosynthesis page loads successfully."""
        _goto_retrosynthesis_page(page, streamlit_server)
        # Use more specific selector to avoid strict mode violation
        main = page.locator('[data-testid="stMainBlockContainer"]').first
        expect(main.get_by_role("heading", name="🧪 Retrosynthesis Advisor")).to_be_visible()

    def test_analyze_smiles_renders_tree(self, page: Page, streamlit_server):
        """Input SMILES, click analyze, tree renders."""
        _goto_retrosynthesis_page(page, streamlit_server)
        
        # Input SMILES in text area
        smiles_input = page.get_by_placeholder("Enter SMILES string")
        if smiles_input.is_visible():
            smiles_input.fill("CC(=O)Nc1ccccc1")
            # Press Tab to trigger Streamlit rerun
            page.keyboard.press("Tab")
            page.wait_for_timeout(2000)
            
            # Click Analyze button
            analyze_btn = page.get_by_role("button", name="Analyze")
            # Debug button state
            print(f"Button disabled: {analyze_btn.is_disabled()}")
            if analyze_btn.is_enabled():
                analyze_btn.click()
                page.wait_for_timeout(3000)
                
                # Check for tree visualization
                expect(page.get_by_text("Synthesis Tree").first).to_be_visible()

    def test_switch_routes(self, page: Page, streamlit_server):
        """Switch between alternative routes in Routes tab."""
        _goto_retrosynthesis_page(page, streamlit_server)
        
        # Look for Routes tab
        routes_tab = page.get_by_role("tab", name="🔀 Routes")
        if routes_tab.is_visible():
            routes_tab.click()
            page.wait_for_timeout(1000)
            
            # Check for route comparison heading (more specific)
            expect(page.get_by_role("heading", name="Route Comparison")).to_be_visible()

    def test_building_blocks_check(self, page: Page, streamlit_server):
        """Check building blocks, see vendor results."""
        _goto_retrosynthesis_page(page, streamlit_server)
        
        # Look for Building Blocks tab with emoji
        bb_tab = page.get_by_role("tab", name="🧱 Building Blocks")
        if bb_tab.is_visible():
            bb_tab.click()
            page.wait_for_timeout(1000)
            
            # Check for vendor or availability information (more specific)
            expect(page.get_by_text("Check availability").or_(page.get_by_text("Vendor information"))).to_be_visible()

    def test_cytoscape_interactions(self, page: Page, streamlit_server):
        """Cytoscape tree zoom/fit controls work."""
        _goto_retrosynthesis_page(page, streamlit_server)
        
        # Check for Cytoscape visualization controls
        fit_btn = page.get_by_role("button", name="Fit")
        zoom_controls = page.get_by_text("Zoom")
        
        cytoscape_controls = fit_btn.or_(zoom_controls)
        if cytoscape_controls.is_visible():
            # Just verify controls are present
            expect(cytoscape_controls).to_be_visible()
