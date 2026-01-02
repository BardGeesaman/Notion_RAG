"""E2E tests for Alignments dashboard."""

import os
import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def alignments_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to Alignments page."""
    os.environ["DISABLE_AUTH"] = "true"
    page.goto(f"{streamlit_server}/?page=Alignments", wait_until="domcontentloaded")
    page.wait_for_timeout(2000)
    return page


class TestAlignmentsE2E:
    """E2E tests for Alignments dashboard page."""

    def test_alignments_page_loads(self, alignments_page: Page):
        """Page loads with all 4 tabs visible."""
        expect(alignments_page.get_by_text("Alignment Files")).to_be_visible()
        expect(alignments_page.get_by_role("tab", name="Browse")).to_be_visible()
        expect(alignments_page.get_by_role("tab", name="Upload")).to_be_visible()
        expect(alignments_page.get_by_role("tab", name="View")).to_be_visible()
        expect(alignments_page.get_by_role("tab", name="Stats")).to_be_visible()

    def test_alignments_upload_tab_has_form(self, alignments_page: Page):
        """Upload tab has file uploaders."""
        upload_tab = alignments_page.get_by_role("tab", name="Upload")
        upload_tab.click()
        alignments_page.wait_for_timeout(500)
        
        expect(alignments_page.get_by_text("Upload Alignment File")).to_be_visible()
        expect(alignments_page.get_by_text("Alignment File (BAM/CRAM)")).to_be_visible()

    def test_alignments_view_tab_has_region_input(self, alignments_page: Page):
        """View tab has region input field."""
        view_tab = alignments_page.get_by_role("tab", name="View")
        view_tab.click()
        alignments_page.wait_for_timeout(500)
        
        expect(alignments_page.get_by_text("View Alignment Reads")).to_be_visible()
        expect(alignments_page.get_by_text("Genomic Region")).to_be_visible()

    def test_alignments_stats_tab_shows_info(self, alignments_page: Page):
        """Stats tab shows guidance when no alignment selected."""
        stats_tab = alignments_page.get_by_role("tab", name="Stats")
        stats_tab.click()
        alignments_page.wait_for_timeout(500)
        
        expect(alignments_page.get_by_text("Alignment Statistics")).to_be_visible()
        expect(alignments_page.get_by_text("Select an alignment")).to_be_visible()
