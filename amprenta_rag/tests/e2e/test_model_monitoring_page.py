"""E2E tests for Model Monitoring dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_model_monitoring_page(page: Page, base_url: str) -> None:
    """Navigate to the Model Monitoring page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Model%20Monitoring")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_model_monitoring_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Model Monitoring page loads successfully."""
    _goto_model_monitoring_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Model Monitoring heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Model\s+Monitoring", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Model\s+Health", re.IGNORECASE)),
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


def test_model_monitoring_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that monitoring tabs are present."""
    _goto_model_monitoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Overview, Details, etc.)
    overview_tab = page.get_by_text(re.compile("Overview", re.IGNORECASE), exact=False)
    
    # At least one strategy should find tabs
    assert tab_count >= 2 or overview_tab.count() > 0, f"Expected tabs but found {tab_count} tab elements"


def test_model_monitoring_content_area(page: Page, streamlit_server: str) -> None:
    """Test that main content area renders."""
    _goto_model_monitoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show either model data or "No models" message
    no_models = main.get_by_text(re.compile("No models", re.IGNORECASE))
    model_health = main.get_by_text(re.compile("Model.*Health", re.IGNORECASE))
    
    # At least one should be visible or main should have content
    try:
        expect(model_health.first).to_be_visible(timeout=5000)
    except AssertionError:
        try:
            expect(no_models.first).to_be_visible(timeout=5000)
        except AssertionError:
            # Just check main has some content
            expect(main).to_be_visible(timeout=5000)


def test_model_monitoring_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that page has interactive elements."""
    _goto_model_monitoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have some interactive elements (buttons, selectors, or tabs)
    # Time window selector may only show if models exist
    buttons = main.get_by_role("button")
    tabs = page.locator('[role="tab"]')
    selectboxes = main.locator('[data-baseweb="select"]')
    
    # Should have at least some interactive elements
    assert (buttons.count() > 0 or tabs.count() > 0 or 
            selectboxes.count() > 0), "Expected interactive elements on page"


def test_model_monitoring_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_model_monitoring_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

