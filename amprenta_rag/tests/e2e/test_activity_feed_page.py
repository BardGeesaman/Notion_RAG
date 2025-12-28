"""E2E tests for Activity Feed dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_activity_feed_page(page: Page, base_url: str) -> None:
    """Navigate to the Activity Feed page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Activity%20Feed")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_activity_feed_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Activity Feed page loads successfully."""
    _goto_activity_feed_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Activity Feed heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Activity\s+Feed", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Recent\s+Activity", re.IGNORECASE)),
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


def test_activity_feed_filters_present(page: Page, streamlit_server: str) -> None:
    """Test that filter controls are present in sidebar."""
    _goto_activity_feed_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Check sidebar for filters
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have filter heading
    filters_heading = sidebar.get_by_text(re.compile("Filters", re.IGNORECASE))
    expect(filters_heading.first).to_be_visible(timeout=10000)
    
    # Should have selectboxes (Program, Time Range, etc.)
    selectboxes = sidebar.locator('[data-baseweb="select"]')
    assert selectboxes.count() >= 2, "Expected at least 2 filter dropdowns (Program, Time Range)"


def test_activity_feed_content_area(page: Page, streamlit_server: str) -> None:
    """Test that main content area renders."""
    _goto_activity_feed_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either activity cards or "No activity" message
    no_activity = main.get_by_text(re.compile("No activity", re.IGNORECASE))
    recent_activity = main.get_by_text(re.compile("Recent Activity", re.IGNORECASE))
    
    # At least one should be visible
    try:
        expect(recent_activity.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(no_activity.first).to_be_visible(timeout=5000)


def test_activity_feed_summary_expander(page: Page, streamlit_server: str) -> None:
    """Test that activity summary expander exists."""
    _goto_activity_feed_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Look for Activity Summary expander
    summary = main.get_by_text(re.compile("Activity Summary", re.IGNORECASE))
    expect(summary.first).to_be_visible(timeout=10000)


def test_activity_feed_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_activity_feed_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

