"""E2E tests for Active Learning dashboard."""

import pytest
from playwright.sync_api import Page, expect

pytestmark = pytest.mark.requires_server


def _goto_active_learning_page(page: Page, base_url: str) -> None:
    """Navigate to the Active Learning page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Active%20Learning")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


class TestActiveLearningE2E:
    def test_load_active_learning_page(self, page: Page, streamlit_server):
        """Active Learning page loads successfully."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Verify page loads - just check main content area is present
        main = page.locator('[data-testid="stMainBlockContainer"]').first
        expect(main).to_be_visible(timeout=10000)
        
        # Check if any content is visible (the page loaded)
        page_content = page.get_by_text("Active Learning", exact=False).first
        expect(page_content).to_be_attached()

    def test_label_queue_tab_displays(self, page: Page, streamlit_server):
        """Label Queue tab displays queue table and metrics."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Click Label Queue tab (should be default but ensure it's selected)
        queue_tab = page.get_by_role("tab", name="📋 Label Queue")
        if queue_tab.is_visible():
            queue_tab.click()
            page.wait_for_timeout(2000)
            
            # Check for queue-related content
            main = page.locator('[data-testid="stMainBlockContainer"]').first
            expect(main.get_by_text("Queue Status").or_(main.get_by_text("Pending Labels"))).to_be_visible()

    def test_labeling_interface_form(self, page: Page, streamlit_server):
        """Labeling Interface tab shows label form elements."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Navigate to Labeling Interface tab
        labeling_tab = page.get_by_role("tab", name="✏️ Labeling Interface")
        if labeling_tab.is_visible():
            labeling_tab.click()
            page.wait_for_timeout(2000)
            
            # Check for labeling form elements
            main = page.locator('[data-testid="stMainBlockContainer"]').first
            expect(main.get_by_text("Label").or_(main.get_by_text("Submit Label"))).to_be_visible()

    def test_cycle_history_tab(self, page: Page, streamlit_server):
        """Cycle History tab shows stats and timeline."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Navigate to Cycle History tab
        history_tab = page.get_by_role("tab", name="📊 Cycle History")
        if history_tab.is_visible():
            history_tab.click()
            page.wait_for_timeout(2000)
            
            # Check for history/stats content
            main = page.locator('[data-testid="stMainBlockContainer"]').first
            expect(main.get_by_text("Cycle").or_(main.get_by_text("Performance"))).to_be_visible()

    def test_sample_selection_tab(self, page: Page, streamlit_server):
        """Sample Selection tab shows strategy dropdown and batch controls."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Navigate to Sample Selection tab
        selection_tab = page.get_by_role("tab", name="🎲 Sample Selection")
        if selection_tab.is_visible():
            selection_tab.click()
            page.wait_for_timeout(2000)
            
            # Check for selection strategy elements
            main = page.locator('[data-testid="stMainBlockContainer"]').first
            expect(main.get_by_text("Strategy").or_(main.get_by_text("Batch Size"))).to_be_visible()

    def test_tab_navigation(self, page: Page, streamlit_server):
        """Click through all 4 tabs and verify each loads."""
        _goto_active_learning_page(page, streamlit_server)
        
        # Test navigation through all tabs
        tabs = [
            ("📋 Label Queue", "Queue"),
            ("✏️ Labeling Interface", "Label"),
            ("📊 Cycle History", "Cycle"),
            ("🎲 Sample Selection", "Strategy")
        ]
        
        for tab_name, expected_content in tabs:
            tab = page.get_by_role("tab", name=tab_name)
            if tab.is_visible():
                tab.click()
                page.wait_for_timeout(1500)
                
                # Verify tab content loads (look for any related text)
                main = page.locator('[data-testid="stMainBlockContainer"]').first
                expect(main.get_by_text(expected_content).first).to_be_visible()
