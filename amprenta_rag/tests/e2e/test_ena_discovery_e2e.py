"""E2E tests for ENA Discovery dashboard."""

import os
import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def ena_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to ENA Discovery page."""
    os.environ["DISABLE_AUTH"] = "true"
    page.goto(f"{streamlit_server}/?page=ENA+Discovery", wait_until="domcontentloaded")
    page.wait_for_timeout(2000)
    return page


class TestENADiscoveryE2E:
    """E2E tests for ENA Discovery page."""

    def test_ena_page_loads(self, ena_page: Page):
        """Page loads with 3 tabs visible."""
        expect(ena_page.get_by_text("ENA Discovery")).to_be_visible()
        expect(ena_page.get_by_role("tab", name="Search")).to_be_visible()
        expect(ena_page.get_by_role("tab", name="Results")).to_be_visible()
        expect(ena_page.get_by_role("tab", name="Ingest")).to_be_visible()

    def test_ena_search_has_controls(self, ena_page: Page):
        """Search tab has keyword input and search button."""
        search_tab = ena_page.get_by_role("tab", name="Search")
        search_tab.click()
        ena_page.wait_for_timeout(500)
        
        expect(ena_page.get_by_placeholder("Homo sapiens")).to_be_visible()
        expect(ena_page.get_by_role("button", name="Search ENA")).to_be_visible()

    def test_ena_results_tab_shows_info(self, ena_page: Page):
        """Results tab shows guidance when empty."""
        results_tab = ena_page.get_by_role("tab", name="Results")
        results_tab.click()
        ena_page.wait_for_timeout(500)
        
        expect(ena_page.get_by_text("No results yet")).to_be_visible()

    def test_ena_ingest_tab_shows_info(self, ena_page: Page):
        """Ingest tab shows guidance when no studies selected."""
        ingest_tab = ena_page.get_by_role("tab", name="Ingest")
        ingest_tab.click()
        ena_page.wait_for_timeout(500)
        
        expect(ena_page.get_by_text("No studies selected")).to_be_visible()

    def test_ena_organism_filter_exists(self, ena_page: Page):
        """Organism filter is available."""
        search_tab = ena_page.get_by_role("tab", name="Search")
        search_tab.click()
        ena_page.wait_for_timeout(500)
        
        expect(ena_page.get_by_text("Organism Filter")).to_be_visible()


@pytest.fixture
def variants_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to Variants page."""
    os.environ["DISABLE_AUTH"] = "true"
    page.goto(f"{streamlit_server}/?page=Variant+Tracking", wait_until="domcontentloaded")
    page.wait_for_timeout(2000)
    return page


class TestVariantsE2E:
    """E2E tests for Variants page with VCF tab."""

    def test_variants_page_loads(self, variants_page: Page):
        """Page loads with 3 tabs visible."""
        expect(variants_page.get_by_text("Variant Tracking")).to_be_visible()
        expect(variants_page.get_by_role("tab", name="Browse")).to_be_visible()
        expect(variants_page.get_by_role("tab", name="Add Variant")).to_be_visible()
        expect(variants_page.get_by_role("tab", name="Upload VCF")).to_be_visible()

    def test_variants_vcf_tab_has_uploader(self, variants_page: Page):
        """VCF tab has file uploader."""
        vcf_tab = variants_page.get_by_role("tab", name="Upload VCF")
        vcf_tab.click()
        variants_page.wait_for_timeout(500)
        
        expect(variants_page.get_by_text("Upload VCF File")).to_be_visible()
        expect(variants_page.get_by_text("Choose VCF file")).to_be_visible()

    def test_variants_browse_has_filters(self, variants_page: Page):
        """Browse tab has gene and organism filters."""
        browse_tab = variants_page.get_by_role("tab", name="Browse")
        browse_tab.click()
        variants_page.wait_for_timeout(500)
        
        expect(variants_page.get_by_placeholder("Filter by gene name")).to_be_visible()
