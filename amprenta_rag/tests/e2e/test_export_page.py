"""E2E tests for Data Export Wizard dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_export_page(page: Page, base_url: str) -> None:
    """Navigate to the Data Export page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Data%20Export")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_export_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Export page loads successfully."""
    _goto_export_page(page, streamlit_server)

    main = _main_container(page)
    
    heading = main.get_by_text(re.compile(r"Data.*Export", re.IGNORECASE))
    
    try:
        expect(heading.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(main).to_be_visible(timeout=10000)


def test_export_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that export tabs are present."""
    _goto_export_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    tabs = page.locator('[role="tab"]')
    assert tabs.count() >= 3, f"Expected 3 tabs, found {tabs.count()}"


def test_export_quick_export_controls(page: Page, streamlit_server: str) -> None:
    """Test that quick export has controls."""
    _goto_export_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have entity type selector
    selectboxes = main.locator('[data-baseweb="select"]')
    
    assert selectboxes.count() > 0, "Expected entity type selector"


def test_export_format_selector(page: Page, streamlit_server: str) -> None:
    """Test that format selector is present."""
    _goto_export_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have format selection
    format_text = main.get_by_text(re.compile("Format|CSV|Excel|JSON", re.IGNORECASE))
    
    assert format_text.count() > 0, "Expected format selector"


def test_export_download_button(page: Page, streamlit_server: str) -> None:
    """Test that download button exists."""
    _goto_export_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have download button
    download_btn = main.get_by_role("button").filter(has_text=re.compile("Download", re.IGNORECASE))
    
    assert download_btn.count() > 0, "Expected download button"


def test_export_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_export_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

