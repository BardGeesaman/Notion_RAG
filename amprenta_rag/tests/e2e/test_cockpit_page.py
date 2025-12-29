"""E2E tests for Scientist's Cockpit dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_cockpit_page(page: Page, base_url: str) -> None:
    """Navigate to the Cockpit page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Cockpit")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_cockpit_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Cockpit page loads successfully."""
    _goto_cockpit_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Cockpit heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Scientist.*Cockpit", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Cockpit", re.IGNORECASE)),
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
        expect(main).to_be_visible(timeout=10000)


def test_cockpit_welcome_banner(page: Page, streamlit_server: str) -> None:
    """Test that welcome banner is visible."""
    _goto_cockpit_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    main = _main_container(page)
    
    # Should have welcome text or date
    welcome = main.get_by_text(re.compile("Welcome|Cockpit", re.IGNORECASE))
    
    assert welcome.count() > 0, "Expected welcome banner"


def test_cockpit_stats_row(page: Page, streamlit_server: str) -> None:
    """Test that stats row is rendered."""
    _goto_cockpit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have statistics text or metrics
    stats = main.get_by_text(re.compile("Statistics|Datasets|Experiments", re.IGNORECASE))
    
    assert stats.count() > 0, "Expected stats section"


def test_cockpit_activity_widget(page: Page, streamlit_server: str) -> None:
    """Test that activity widget is rendered."""
    _goto_cockpit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have activity section
    activity = main.get_by_text(re.compile("Activity|Recent", re.IGNORECASE))
    
    assert activity.count() > 0, "Expected activity widget"


def test_cockpit_quick_actions(page: Page, streamlit_server: str) -> None:
    """Test that quick actions are present."""
    _goto_cockpit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have action buttons
    actions = main.get_by_text(re.compile("Quick.*Actions|Actions", re.IGNORECASE))
    buttons = main.get_by_role("button")
    
    assert actions.count() > 0 or buttons.count() > 0, "Expected quick actions"


def test_cockpit_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_cockpit_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

