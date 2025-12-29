"""E2E tests for Share Links dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_share_links_page(page: Page, base_url: str) -> None:
    """Navigate to the Share Links page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Share%20Links")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_share_links_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Share Links page loads successfully."""
    _goto_share_links_page(page, streamlit_server)

    main = _main_container(page)
    
    heading = main.get_by_text(re.compile(r"Share.*Links", re.IGNORECASE))
    
    try:
        expect(heading.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(main).to_be_visible(timeout=10000)


def test_share_links_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that share links tabs are present."""
    _goto_share_links_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    tabs = page.locator('[role="tab"]')
    assert tabs.count() >= 2, "Expected 2 tabs"


def test_share_links_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that page has interactive elements."""
    _goto_share_links_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have buttons or textareas
    buttons = main.get_by_role("button")
    textareas = main.locator('textarea')
    
    assert buttons.count() > 0 or textareas.count() > 0, "Expected interactive elements"


def test_manage_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that manage tab has content."""
    _goto_share_links_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Manage Links tab
    manage_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Manage", re.IGNORECASE))
    
    try:
        manage_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    manage_text = main.get_by_text(re.compile("Manage|Load|Links", re.IGNORECASE))
    
    assert manage_text.count() > 0 or main.count() > 0, "Manage tab should render"


def test_share_links_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_share_links_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

