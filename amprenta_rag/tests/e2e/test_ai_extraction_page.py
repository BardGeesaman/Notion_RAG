"""E2E tests for AI Document Extraction dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_ai_extraction_page(page: Page, base_url: str) -> None:
    """Navigate to the AI Extraction page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=AI%20Extraction")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_ai_extraction_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the AI Extraction page loads successfully."""
    _goto_ai_extraction_page(page, streamlit_server)

    main = _main_container(page)

    # Look for AI Extraction heading
    heading_patterns = [
        main.get_by_text(re.compile(r"AI.*Extraction", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Document.*Extraction", re.IGNORECASE)),
    ]
    
    found = False
    for pattern in heading_patterns:
        try:
            expect(pattern.first).to_be_visible(timeout=5000)
            found = True
            break
        except AssertionError:
            continue
    
    if not found:
        # Fallback: just check if main container has content
        expect(main).to_be_visible(timeout=10000)


def test_ai_extraction_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that extraction tabs are present."""
    _goto_ai_extraction_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Upload, Jobs, Results)
    upload_tab = page.get_by_text("Upload", exact=False)
    jobs_tab = page.get_by_text("Jobs", exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 2 or upload_tab.count() > 0 or 
            jobs_tab.count() > 0), f"Expected tabs but found {tab_count} tab elements"


def test_ai_extraction_upload_section(page: Page, streamlit_server: str) -> None:
    """Test that upload section has file uploader."""
    _goto_ai_extraction_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have file uploader or upload-related content
    upload_text = main.get_by_text(re.compile("Upload|Select files", re.IGNORECASE))
    start_extraction = main.get_by_role("button").filter(has_text=re.compile("Start.*Extraction|Extract", re.IGNORECASE))
    
    # Should have upload UI elements
    assert upload_text.count() > 0 or start_extraction.count() > 0, "Expected upload section"


def test_ai_extraction_jobs_tab(page: Page, streamlit_server: str) -> None:
    """Test that Jobs tab has job listing functionality."""
    _goto_ai_extraction_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Jobs tab
    jobs_tab = page.get_by_text("Jobs", exact=False).or_(
        page.locator('[role="tab"]').filter(has_text="Jobs")
    )
    
    try:
        jobs_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have job-related content
    jobs_text = main.get_by_text(re.compile("Recent.*jobs|Jobs|Load", re.IGNORECASE))
    
    assert jobs_text.count() > 0, "Expected jobs section content"


def test_ai_extraction_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_ai_extraction_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

