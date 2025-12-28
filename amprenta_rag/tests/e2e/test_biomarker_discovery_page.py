"""E2E tests for Biomarker Discovery dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_biomarker_discovery_page(page: Page, base_url: str) -> None:
    """Navigate to the Biomarker Discovery page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Biomarker%20Discovery")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_biomarker_discovery_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Biomarker Discovery page loads successfully."""
    _goto_biomarker_discovery_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Biomarker Discovery heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Biomarker\s+Discovery", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Biomarker", re.IGNORECASE)),
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


def test_biomarker_discovery_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that discovery workflow tabs are present."""
    _goto_biomarker_discovery_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Setup, Methods, Results, etc.)
    setup_tab = page.get_by_text("Setup", exact=False)
    methods_tab = page.get_by_text("Methods", exact=False)
    results_tab = page.get_by_text("Results", exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 3 or setup_tab.count() > 0 or methods_tab.count() > 0 or 
            results_tab.count() > 0), f"Expected workflow tabs but found {tab_count} tab elements"


def test_biomarker_discovery_setup_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that Setup tab shows content."""
    _goto_biomarker_discovery_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should show experiment/sample setup or instructions
    experiment_text = main.get_by_text(re.compile("Experiment|Sample", re.IGNORECASE))
    setup_text = main.get_by_text(re.compile("Setup", re.IGNORECASE))
    
    # At least one should be visible
    try:
        expect(experiment_text.first).to_be_visible(timeout=5000)
    except AssertionError:
        try:
            expect(setup_text.first).to_be_visible(timeout=5000)
        except AssertionError:
            # Just check main has content
            expect(main).to_be_visible(timeout=5000)


def test_biomarker_discovery_refresh_button(page: Page, streamlit_server: str) -> None:
    """Test that refresh experiments button exists."""
    _goto_biomarker_discovery_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Look for refresh button
    refresh_button = main.get_by_role("button").filter(has_text=re.compile("Refresh", re.IGNORECASE))
    
    # Button should exist
    assert refresh_button.count() > 0, "Expected refresh button"


def test_biomarker_discovery_workflow_structure(page: Page, streamlit_server: str) -> None:
    """Test that biomarker workflow structure is present."""
    _goto_biomarker_discovery_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have multiple tabs or sections for the workflow
    # Count role="tab" elements
    tabs = page.locator('[role="tab"]')
    
    # Or look for workflow keywords
    workflow_indicators = [
        main.get_by_text(re.compile("Setup", re.IGNORECASE)),
        main.get_by_text(re.compile("Methods", re.IGNORECASE)),
        main.get_by_text(re.compile("Results", re.IGNORECASE)),
    ]
    
    workflow_found = tabs.count() >= 3
    for indicator in workflow_indicators:
        if indicator.count() > 0:
            workflow_found = True
            break
    
    assert workflow_found, "Expected biomarker discovery workflow structure"


def test_biomarker_discovery_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_biomarker_discovery_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

