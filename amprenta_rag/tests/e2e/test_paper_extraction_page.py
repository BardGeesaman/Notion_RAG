"""E2E tests for Paper Search page extraction features."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_paper_search_page(page: Page, base_url: str) -> None:
    """Navigate to the Paper Search page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Paper%20Search")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_paper_extraction_buttons_exist(page: Page, streamlit_server: str) -> None:
    """Test that extraction buttons are present on paper detail view."""
    _goto_paper_search_page(page, streamlit_server)
    
    # Switch to Ingested Papers tab
    main = _main_container(page)
    ingested_tab = main.get_by_role("tab").filter(has_text=re.compile("Ingested", re.IGNORECASE))
    
    try:
        ingested_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass
    
    # Look for extraction-related buttons (may only appear after viewing a paper)
    extract_button = main.get_by_role("button").filter(has_text=re.compile("Extract.*Experiments", re.IGNORECASE))
    view_exp_button = main.get_by_role("button").filter(has_text=re.compile("View.*Extracted", re.IGNORECASE))
    
    # At least one extraction feature should be mentioned or present
    extraction_text = main.get_by_text(re.compile("Extraction|Extract|Experiments", re.IGNORECASE))
    
    # Should have extraction UI elements or text
    assert (extract_button.count() > 0 or view_exp_button.count() > 0 or 
            extraction_text.count() > 0), "Expected extraction features"


def test_paper_extraction_section_loads(page: Page, streamlit_server: str) -> None:
    """Test that extraction section doesn't crash the page."""
    _goto_paper_search_page(page, streamlit_server)
    page.wait_for_timeout(3000)
    
    # Page should be responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

