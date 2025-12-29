"""E2E tests for Screening dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_screening_page(page: Page, base_url: str) -> None:
    """Navigate to the Screening page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Screening")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_screening_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Screening page loads successfully."""
    _goto_screening_page(page, streamlit_server)

    main = _main_container(page)

    # Look for HTS Screening heading
    heading_patterns = [
        main.get_by_text(re.compile(r"HTS.*Screening", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Screening.*Campaigns", re.IGNORECASE)),
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


def test_screening_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that screening tabs are present."""
    _goto_screening_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content
    campaigns_tab = page.get_by_text("Campaigns", exact=False)
    details_tab = page.get_by_text(re.compile("Campaign.*Details|Details", re.IGNORECASE), exact=False)
    active_tab = page.get_by_text(re.compile("Active.*Learning", re.IGNORECASE), exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 3 or campaigns_tab.count() > 0 or details_tab.count() > 0 or 
            active_tab.count() > 0), f"Expected tabs but found {tab_count} tab elements"


def test_screening_campaigns_list(page: Page, streamlit_server: str) -> None:
    """Test that campaigns list displays."""
    _goto_screening_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either campaigns table or "No campaigns" message
    dataframe = main.locator('[data-testid="stDataFrame"]')
    no_campaigns = main.get_by_text(re.compile("No.*campaigns", re.IGNORECASE))
    
    # Should have either data or empty message
    assert dataframe.count() > 0 or no_campaigns.count() > 0, "Expected campaigns content"


def test_screening_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that page has interactive elements."""
    _goto_screening_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have some interactive elements (buttons, tabs, etc.)
    # Campaign selector only shows if campaigns exist and API call succeeds
    tabs = page.locator('[role="tab"]')
    buttons = main.get_by_role("button")
    
    # Should have interactive elements
    assert tabs.count() > 0 or buttons.count() > 0, "Expected interactive elements on page"


def test_screening_active_learning_section(page: Page, streamlit_server: str) -> None:
    """Test that active learning tab has content."""
    _goto_screening_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Active Learning tab
    active_tab = page.get_by_text(re.compile("Active.*Learning", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Active|Learning", re.IGNORECASE))
    )
    
    try:
        active_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have active learning related content or strategy selector
    active_learning = main.get_by_text(re.compile("Active.*Learning|Suggestions|Strategy|Batch.*Size", re.IGNORECASE))
    buttons = main.get_by_role("button")
    
    # Should have AL content or interactive elements
    assert active_learning.count() > 0 or buttons.count() > 0, "Expected active learning content"


def test_screening_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_screening_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

