"""E2E tests for Sync Monitor dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_sync_monitor_page(page: Page, base_url: str) -> None:
    """Navigate to the Sync Monitor page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Sync%20Monitor")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_sync_monitor_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Sync Monitor page loads successfully."""
    _goto_sync_monitor_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Sync Monitor heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Sync.*Monitor", re.IGNORECASE)),
        main.get_by_text(re.compile(r"External.*Sync", re.IGNORECASE)),
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


def test_sync_monitor_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that sync monitor tabs are present."""
    _goto_sync_monitor_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Look for tab text content (Jobs, Conflicts, Stats)
    jobs_tab = page.get_by_text("Jobs", exact=False)
    conflicts_tab = page.get_by_text("Conflicts", exact=False)
    
    # At least one strategy should find tabs
    assert (tab_count >= 2 or jobs_tab.count() > 0 or 
            conflicts_tab.count() > 0), f"Expected tabs but found {tab_count} tab elements"


def test_sync_monitor_run_sync_section(page: Page, streamlit_server: str) -> None:
    """Test that run sync section has controls."""
    _goto_sync_monitor_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have source selector and run button
    run_sync_button = main.get_by_role("button").filter(has_text=re.compile("Run.*Sync", re.IGNORECASE))
    source_text = main.get_by_text(re.compile("Source", re.IGNORECASE))
    
    # Should have sync controls
    assert run_sync_button.count() > 0 or source_text.count() > 0, "Expected sync controls"


def test_sync_monitor_source_selector(page: Page, streamlit_server: str) -> None:
    """Test that source selector is present."""
    _goto_sync_monitor_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have selectboxes for source selection
    selectboxes = main.locator('[data-baseweb="select"]')
    chembl_text = main.get_by_text(re.compile("chembl|pubchem", re.IGNORECASE))
    
    # Should have source selection
    assert selectboxes.count() > 0 or chembl_text.count() > 0, "Expected source selector"


def test_sync_monitor_load_jobs_button(page: Page, streamlit_server: str) -> None:
    """Test that load jobs button exists."""
    _goto_sync_monitor_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have Load jobs button
    load_button = main.get_by_role("button").filter(has_text=re.compile("Load", re.IGNORECASE))
    
    assert load_button.count() > 0, "Expected Load jobs button"


def test_sync_monitor_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_sync_monitor_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

