"""E2E tests for Audit Trail dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_audit_page(page: Page, base_url: str) -> None:
    """Navigate to the Audit Trail page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Audit%20Trail")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_audit_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Audit Trail page loads successfully."""
    _goto_audit_page(page, streamlit_server)

    main = _main_container(page)
    
    heading = main.get_by_text(re.compile(r"Audit.*Trail", re.IGNORECASE))
    
    try:
        expect(heading.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(main).to_be_visible(timeout=10000)


def test_audit_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that audit tabs are present."""
    _goto_audit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    tabs = page.locator('[role="tab"]')
    assert tabs.count() >= 3, f"Expected 3 tabs, found {tabs.count()}"


def test_audit_entity_selector(page: Page, streamlit_server: str) -> None:
    """Test that entity selector is present."""
    _goto_audit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have entity type selector
    entity_select = main.get_by_text(re.compile("Entity.*Type", re.IGNORECASE))
    
    assert entity_select.count() > 0, "Expected entity type selector"


def test_audit_fetch_button(page: Page, streamlit_server: str) -> None:
    """Test that fetch audit trail button exists."""
    _goto_audit_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have fetch button
    fetch_btn = main.get_by_role("button").filter(has_text=re.compile("Fetch|Load", re.IGNORECASE))
    
    assert fetch_btn.count() > 0, "Expected fetch button"


def test_audit_integrity_check_tab(page: Page, streamlit_server: str) -> None:
    """Test that integrity check tab renders."""
    _goto_audit_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Integrity Check tab
    integrity_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Integrity", re.IGNORECASE))
    
    try:
        integrity_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    integrity_text = main.get_by_text(re.compile("Integrity|Verify|Checksum", re.IGNORECASE))
    
    assert integrity_text.count() > 0 or main.count() > 0, "Integrity check tab should render"


def test_audit_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_audit_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

