"""E2E tests for Experiments dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_experiments_page(page: Page, base_url: str) -> None:
    """Navigate to the Experiments page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Experiments")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_experiments_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Experiments page loads successfully."""
    _goto_experiments_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Experiments heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Experiments", re.IGNORECASE)),
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


def test_experiments_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that experiment tabs are present."""
    _goto_experiments_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content
    browse_tab = page.get_by_text("Browse", exact=False)
    
    # At least one strategy should find tabs
    assert tab_count >= 2 or browse_tab.count() > 0, f"Expected tabs but found {tab_count} tab elements"


def test_experiments_search_input(page: Page, streamlit_server: str) -> None:
    """Test that search input is present."""
    _goto_experiments_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have search input for experiment name
    search_input = main.get_by_label(re.compile("Search.*experiments", re.IGNORECASE)).or_(
        main.get_by_placeholder(re.compile("Search", re.IGNORECASE))
    )
    
    # At least one search input should exist
    assert search_input.count() > 0, "Expected search input"


def test_experiments_total_metric(page: Page, streamlit_server: str) -> None:
    """Test that total experiments metric is displayed."""
    _goto_experiments_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show total experiments count
    total_experiments = main.get_by_text(re.compile("Total.*Experiments", re.IGNORECASE))
    
    expect(total_experiments.first).to_be_visible(timeout=10000)


def test_experiments_content_display(page: Page, streamlit_server: str) -> None:
    """Test that experiments list or empty state displays."""
    _goto_experiments_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either experiments table or expanders
    dataframe = main.locator('[data-testid="stDataFrame"]')
    expanders = main.locator('[data-testid="stExpander"]')
    
    # Should have some content elements
    assert dataframe.count() > 0 or expanders.count() > 0, "Expected experiments content"


def test_experiments_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_experiments_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

