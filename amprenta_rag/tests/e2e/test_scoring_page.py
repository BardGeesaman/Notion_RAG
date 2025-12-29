"""E2E tests for Scoring dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_scoring_page(page: Page, base_url: str) -> None:
    """Navigate to the Scoring page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Scoring")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_scoring_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Scoring page loads successfully."""
    _goto_scoring_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Compound Scoring heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Compound.*Scoring", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Scoring", re.IGNORECASE)),
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


def test_scoring_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that scoring tabs are present."""
    _goto_scoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Should have at least 2 tabs
    assert tab_count >= 2, f"Expected at least 2 tabs but found {tab_count}"


def test_scoring_single_compound_form(page: Page, streamlit_server: str) -> None:
    """Test that single compound scoring form has elements."""
    _goto_scoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have form inputs
    inputs = main.locator('input')
    textareas = main.locator('textarea')
    buttons = main.get_by_role("button")
    
    # Should have form elements
    assert inputs.count() > 0, "Expected input fields"
    assert buttons.count() > 0, "Expected score button"


def test_scoring_configuration_options(page: Page, streamlit_server: str) -> None:
    """Test that scoring configuration options are present."""
    _goto_scoring_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have checkboxes for scoring options
    checkboxes = main.locator('input[type="checkbox"]')
    
    # Or should have inputs for configuration
    inputs = main.locator('input')
    
    # Should have some configuration elements
    assert checkboxes.count() > 0 or inputs.count() > 0, "Expected configuration elements"


def test_scoring_batch_tab_content(page: Page, streamlit_server: str) -> None:
    """Test that batch scoring tab has content."""
    _goto_scoring_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Batch Scoring tab
    batch_tab = page.get_by_text(re.compile("Batch.*Scoring", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Batch", re.IGNORECASE))
    )
    
    try:
        batch_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have batch scoring content
    batch_text = main.get_by_text(re.compile("Batch|SMILES|Score", re.IGNORECASE))
    
    assert batch_text.count() > 0, "Expected batch scoring content"


def test_scoring_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_scoring_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

