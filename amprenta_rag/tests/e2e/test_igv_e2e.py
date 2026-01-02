"""E2E tests for IGV.js genome browser integration."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def base_url(request):
    """Get base URL from pytest config or default."""
    return getattr(request.config, "base_url", "http://localhost:8501")


@pytest.mark.requires_server
class TestIGVGenomeBrowser:
    """E2E tests for IGV.js genome browser page."""
    
    def test_genome_browser_loads(self, page: Page, base_url: str):
        """Test genome browser page renders."""
        page.goto(f"{base_url}/Genome_Browser", wait_until="domcontentloaded")
        
        # Should show header
        expect(page.get_by_text("IGV.js Genome Browser")).to_be_visible(timeout=10000)
        
        # Should show reference genome selector
        expect(page.get_by_text("Reference genome")).to_be_visible()
    
    def test_alignment_selector_present(self, page: Page, base_url: str):
        """Test alignment file selector is shown."""
        page.goto(f"{base_url}/Genome_Browser", wait_until="domcontentloaded")
        
        # Should show alignment tracks section
        expect(page.get_by_text("Alignment Tracks")).to_be_visible(timeout=10000)
        
        # Should show either multiselect or "no files" message
        alignment_section = page.locator("text=BAM/CRAM").or_(page.locator("text=No indexed alignment"))
        expect(alignment_section).to_be_visible()
    
    def test_cram_reference_input_present(self, page: Page, base_url: str):
        """Test CRAM reference URL input is shown."""
        page.goto(f"{base_url}/Genome_Browser", wait_until="domcontentloaded")
        
        # Should show CRAM reference input
        expect(page.get_by_text("CRAM Reference FASTA URL")).to_be_visible(timeout=10000)
    
    def test_load_browser_button(self, page: Page, base_url: str):
        """Test Load Genome Browser button exists."""
        page.goto(f"{base_url}/Genome_Browser", wait_until="domcontentloaded")
        
        # Should show load button
        expect(page.get_by_role("button", name="Load Genome Browser")).to_be_visible(timeout=10000)
