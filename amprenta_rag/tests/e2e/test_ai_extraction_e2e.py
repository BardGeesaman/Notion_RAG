"""E2E tests for AI Extraction Tools dashboard."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def extraction_page(page: Page, streamlit_server: str) -> Page:
    """Navigate to AI Extraction Tools page."""
    page.goto(f"{streamlit_server}/?page=AI+Extraction")
    page.wait_for_load_state("domcontentloaded")
    expect(page.get_by_role("heading", name="AI Extraction Tools")).to_be_visible(timeout=10000)
    return page


class TestAIExtractionE2E:
    """E2E tests for AI Extraction Tools dashboard."""

    def test_extraction_page_loads(self, extraction_page: Page):
        """Page renders with 3 tabs."""
        expect(extraction_page.get_by_text("OCR")).to_be_visible()
        expect(extraction_page.get_by_text("Web Scraper")).to_be_visible()
        expect(extraction_page.get_by_text("Entity Normalizer")).to_be_visible()

    def test_ocr_tab_has_upload(self, extraction_page: Page):
        """OCR tab has file upload and language selector."""
        extraction_page.get_by_text("OCR").click()
        extraction_page.wait_for_timeout(500)
        
        # Should have file uploader
        expect(extraction_page.get_by_text("Upload file")).to_be_visible()
        
        # Should have language selector
        expect(extraction_page.get_by_text("Language")).to_be_visible()
        
        # Should have extract button (disabled without file)
        expect(extraction_page.get_by_role("button", name="Extract Text")).to_be_visible()

    def test_normalizer_lookup(self, extraction_page: Page):
        """Entity normalizer allows lookup."""
        extraction_page.get_by_text("Entity Normalizer").click()
        extraction_page.wait_for_timeout(500)
        
        # Should have entity type selector
        expect(extraction_page.get_by_text("Entity Type")).to_be_visible()
        
        # Should have name input
        entity_input = extraction_page.get_by_placeholder("e.g., BRCA1, aspirin, diabetes")
        expect(entity_input).to_be_visible()
        
        # Enter a test entity
        entity_input.fill("BRCA1")
        extraction_page.keyboard.press("Tab")
        extraction_page.wait_for_timeout(500)
        
        # Should have lookup button
        expect(extraction_page.get_by_role("button", name="Lookup")).to_be_visible()
