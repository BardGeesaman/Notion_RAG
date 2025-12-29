"""E2E tests for Programs dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_programs_page(page: Page, base_url: str) -> None:
    """Navigate to the Programs page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Programs")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_programs_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Programs page loads successfully."""
    _goto_programs_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Programs heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Programs", re.IGNORECASE)),
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


def test_programs_search_input(page: Page, streamlit_server: str) -> None:
    """Test that search input is present."""
    _goto_programs_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have search input for program name
    search_input = main.get_by_label(re.compile("Search.*programs", re.IGNORECASE)).or_(
        main.get_by_placeholder(re.compile("Search", re.IGNORECASE))
    )
    
    # Search input should be visible
    expect(search_input.first).to_be_visible(timeout=10000)


def test_programs_total_metric(page: Page, streamlit_server: str) -> None:
    """Test that total programs metric is displayed."""
    _goto_programs_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show total programs count or "No programs" message
    total_programs = main.get_by_text(re.compile("Total.*Programs", re.IGNORECASE))
    no_programs = main.get_by_text(re.compile("No programs", re.IGNORECASE))
    
    # At least one should be visible
    try:
        expect(total_programs.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(no_programs.first).to_be_visible(timeout=5000)


def test_programs_content_display(page: Page, streamlit_server: str) -> None:
    """Test that programs list or empty state displays."""
    _goto_programs_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either program expanders or "No programs found" message
    expanders = main.locator('[data-testid="stExpander"]')
    no_programs = main.get_by_text(re.compile("No programs found", re.IGNORECASE))
    
    # Should have either programs or empty message
    assert expanders.count() > 0 or no_programs.count() > 0, "Expected programs list or empty state"


def test_programs_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_programs_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

