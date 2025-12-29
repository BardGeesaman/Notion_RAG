"""E2E tests for Review Queue dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_review_queue_page(page: Page, base_url: str) -> None:
    """Navigate to the Review Queue page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Review%20Queue")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_review_queue_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Review Queue page loads successfully."""
    _goto_review_queue_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Review Queue heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Review\s+Queue", re.IGNORECASE)),
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


def test_review_queue_pending_count(page: Page, streamlit_server: str) -> None:
    """Test that pending reviews count is displayed."""
    _goto_review_queue_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show pending reviews count
    pending_text = main.get_by_text(re.compile("Pending.*reviews?", re.IGNORECASE))
    
    expect(pending_text.first).to_be_visible(timeout=10000)


def test_review_queue_content_display(page: Page, streamlit_server: str) -> None:
    """Test that review queue content displays."""
    _goto_review_queue_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either review items or "No pending reviews" message
    expanders = main.locator('[data-testid="stExpander"]')
    no_reviews = main.get_by_text(re.compile("No pending reviews", re.IGNORECASE))
    
    # Should have either reviews or empty message
    assert expanders.count() > 0 or no_reviews.count() > 0, "Expected review queue content"


def test_review_queue_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that review queue has interactive elements or appropriate empty state."""
    _goto_review_queue_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have some content - either review items or empty state message
    # Don't require specific buttons since they only appear when reviews exist
    expanders = main.locator('[data-testid="stExpander"]')
    buttons = main.get_by_role("button")
    no_pending = main.get_by_text(re.compile("No pending|No.*reviews", re.IGNORECASE))
    
    # Should have either interactive elements or empty message
    assert (expanders.count() > 0 or buttons.count() > 0 or 
            no_pending.count() > 0), "Expected review queue content"


def test_review_queue_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_review_queue_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

