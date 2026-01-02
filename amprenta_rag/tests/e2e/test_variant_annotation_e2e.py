"""End-to-end tests for variant annotation UI."""

import pytest
from playwright.sync_api import Page, expect


class TestVariantAnnotationUI:
    """E2E tests for variant annotation dashboard."""
    
    @pytest.fixture
    def base_url(self):
        """Dashboard base URL."""
        return "http://localhost:8501"
    
    def test_annotations_tab_visible(self, page: Page, base_url: str):
        """Test that Annotations tab is visible on variants page."""
        page.goto(f"{base_url}/Variants")
        page.wait_for_load_state("domcontentloaded")
        
        annotations_tab = page.get_by_text("Annotations", exact=False)
        expect(annotations_tab).to_be_visible()
    
    def test_batch_annotation_section(self, page: Page, base_url: str):
        """Test batch annotation section in Browse tab."""
        page.goto(f"{base_url}/Variants")
        page.wait_for_load_state("domcontentloaded")
        
        page.get_by_text("Browse").click()
        
        expect(page.get_by_text("Batch Annotation")).to_be_visible()
