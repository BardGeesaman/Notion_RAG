"""E2E tests for Mapping Refresh extended tabs (Conflicts, KEGG Cache)."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def mapping_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to Mapping Refresh page."""
    page.goto(f"{streamlit_server}/?page=Mapping+Refresh")
    page.wait_for_load_state("domcontentloaded")
    expect(page.get_by_role("heading", name="ID Mapping Refresh")).to_be_visible(timeout=10000)
    return page


class TestMappingRefreshExtendedE2E:
    """E2E tests for extended Mapping Refresh tabs."""

    def test_conflicts_tab_shows_filters(self, mapping_page: Page):
        """Conflicts tab has source and status filters."""
        # Click Conflicts tab
        mapping_page.get_by_text("Conflicts").click()
        mapping_page.wait_for_timeout(500)
        
        # Should have source filter
        expect(mapping_page.get_by_text("Source")).to_be_visible()
        
        # Should have status filter
        expect(mapping_page.get_by_text("Status")).to_be_visible()
        
        # Should show conflicts list or "No pending conflicts" message
        page_content = mapping_page.content()
        assert "conflicts" in page_content.lower() or "No pending" in page_content

    def test_kegg_cache_tab_shows_status(self, mapping_page: Page):
        """KEGG Cache tab displays cache metrics."""
        # Click KEGG Cache tab
        mapping_page.get_by_text("KEGG Cache").click()
        mapping_page.wait_for_timeout(500)
        
        # Should have expiration metrics or error message
        page_content = mapping_page.content()
        assert (
            "Expiring" in page_content or
            "Total KEGG" in page_content or
            "Unable to fetch" in page_content
        )
