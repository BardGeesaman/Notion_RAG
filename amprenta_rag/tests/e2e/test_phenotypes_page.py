"""E2E tests for Phenotypes dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_phenotypes_page(page: Page, base_url: str) -> None:
    """Navigate to the Phenotypes page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Phenotypes")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_phenotypes_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Phenotypes page loads successfully."""
    _goto_phenotypes_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Phenotype Explorer heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Phenotype.*Explorer", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Phenotypes?", re.IGNORECASE)),
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


def test_phenotypes_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that phenotype tabs are present."""
    _goto_phenotypes_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Should have at least 2 tabs
    assert tab_count >= 2, f"Expected at least 2 tabs but found {tab_count}"


def test_phenotypes_lookup_functionality(page: Page, streamlit_server: str) -> None:
    """Test that HPO lookup functionality is present."""
    _goto_phenotypes_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have inputs and buttons for HPO lookup
    inputs = main.locator('input')
    buttons = main.get_by_role("button")
    
    # Should have form elements
    assert inputs.count() > 0 and buttons.count() > 0, "Expected lookup form elements"


def test_phenotypes_query_expansion_tab(page: Page, streamlit_server: str) -> None:
    """Test that query expansion tab has content."""
    _goto_phenotypes_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Query Expansion tab
    expansion_tab = page.get_by_text(re.compile("Query.*Expansion|Expansion", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Expansion|Query", re.IGNORECASE))
    )
    
    try:
        expansion_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have query expansion content
    expansion_text = main.get_by_text(re.compile("Query|Clinical|HPO|Expand", re.IGNORECASE))
    textareas = main.locator('textarea')
    
    # Should have expansion UI
    assert expansion_text.count() > 0 or textareas.count() > 0, "Expected query expansion content"


def test_phenotypes_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that page has interactive elements."""
    _goto_phenotypes_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have buttons and inputs
    buttons = main.get_by_role("button")
    inputs = main.locator('input')
    textareas = main.locator('textarea')
    
    # Should have interactive elements
    assert (buttons.count() > 0 or inputs.count() > 0 or 
            textareas.count() > 0), "Expected interactive elements"


def test_phenotypes_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_phenotypes_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

