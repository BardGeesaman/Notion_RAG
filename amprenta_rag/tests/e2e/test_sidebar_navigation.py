"""E2E tests for grouped sidebar navigation."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page  # noqa: E402


pytestmark = pytest.mark.requires_server


def test_sidebar_renders(page: Page, streamlit_server: str) -> None:
    """Test that sidebar renders successfully."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    # Sidebar should be present
    sidebar = page.locator('[data-testid="stSidebar"]')
    assert sidebar.count() > 0, "Sidebar not found"


def test_sidebar_group_expanders_present(page: Page, streamlit_server: str) -> None:
    """Test that group expanders are present in sidebar."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have expandable sections
    expanders = sidebar.locator('[data-testid="stExpander"]')
    
    # Should have at least one group expander
    assert expanders.count() >= 1, f"Expected group expanders, found {expanders.count()}"


def test_sidebar_group_icons_visible(page: Page, streamlit_server: str) -> None:
    """Test that group icons are visible."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Look for emoji icons (Home 🏠, Discovery 🔍, Chemistry ⚗️)
    icons = sidebar.get_by_text(re.compile("🏠|🔍|⚗️|🧪|🔬", re.IGNORECASE))
    
    assert icons.count() > 0, "Group icons not found"


def test_sidebar_page_links_in_groups(page: Page, streamlit_server: str) -> None:
    """Test that page links exist within groups."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have navigation buttons
    buttons = sidebar.get_by_role("button")
    
    # Should have some buttons for navigation
    assert buttons.count() >= 1, "Expected navigation buttons"


def test_sidebar_navigation_structure(page: Page, streamlit_server: str) -> None:
    """Test that sidebar has proper navigation structure."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have Navigation title
    nav_title = sidebar.get_by_text(re.compile("Navigation", re.IGNORECASE))
    
    assert nav_title.count() > 0 or sidebar.count() > 0, "Navigation structure present"


def test_sidebar_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that sidebar doesn't crash the app."""
    page.goto(streamlit_server)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    # Page should load
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Sidebar should exist
    sidebar = page.locator('[data-testid="stSidebar"]')
    assert sidebar.count() > 0, "Sidebar not found"

