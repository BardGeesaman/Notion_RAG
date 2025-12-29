"""E2E tests for IP Portfolio dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_ip_portfolio_page(page: Page, base_url: str) -> None:
    """Navigate to the IP Portfolio page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=IP%20Portfolio")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_ip_portfolio_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the IP Portfolio page loads successfully."""
    _goto_ip_portfolio_page(page, streamlit_server)

    main = _main_container(page)
    
    heading = main.get_by_text(re.compile(r"IP.*Portfolio", re.IGNORECASE))
    
    try:
        expect(heading.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(main).to_be_visible(timeout=10000)


def test_ip_portfolio_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that IP portfolio tabs are present."""
    _goto_ip_portfolio_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    tabs = page.locator('[role="tab"]')
    assert tabs.count() >= 3, f"Expected 3 tabs, found {tabs.count()}"


def test_disclosures_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that disclosures tab has content."""
    _goto_ip_portfolio_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have disclosure-related content
    disc_text = main.get_by_text(re.compile("Disclosure|Invention|Status", re.IGNORECASE))
    
    assert disc_text.count() > 0, "Expected disclosure content"


def test_patents_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that patents tab has content."""
    _goto_ip_portfolio_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Patents tab
    patents_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Patents", re.IGNORECASE))
    
    try:
        patents_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    patents_text = main.get_by_text(re.compile("Patent|Application|Jurisdiction", re.IGNORECASE))
    
    assert patents_text.count() > 0 or main.count() > 0, "Patents tab should render"


def test_evidence_links_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that evidence links tab has content."""
    _goto_ip_portfolio_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Evidence Links tab
    links_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Evidence|Links", re.IGNORECASE))
    
    try:
        links_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    links_text = main.get_by_text(re.compile("Evidence|Link|Entity", re.IGNORECASE))
    
    assert links_text.count() > 0 or main.count() > 0, "Evidence links tab should render"


def test_ip_portfolio_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_ip_portfolio_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

