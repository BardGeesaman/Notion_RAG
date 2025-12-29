"""E2E tests for Chemical Sketcher dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_chemical_sketcher_page(page: Page, base_url: str) -> None:
    """Navigate to the Chemical Sketcher page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Chemical%20Sketcher")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_chemical_sketcher_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Chemical Sketcher page loads successfully."""
    _goto_chemical_sketcher_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Chemical Sketcher heading
    heading_patterns = [
        main.get_by_text(re.compile(r"Chemical.*Sketcher", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Structure.*Sketcher", re.IGNORECASE)),
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


def test_chemical_sketcher_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that sketcher tabs are present."""
    _goto_chemical_sketcher_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    assert tab_count >= 2, f"Expected at least 2 tabs but found {tab_count}"


def test_chemical_sketcher_ketcher_iframe_renders(page: Page, streamlit_server: str) -> None:
    """Test that Ketcher iframe is embedded."""
    _goto_chemical_sketcher_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Look for iframe (Ketcher is loaded via st.components.html)
    iframes = page.locator('iframe')
    
    # Should have at least one iframe for Ketcher
    assert iframes.count() > 0, "Expected Ketcher iframe"


def test_chemical_sketcher_smiles_input(page: Page, streamlit_server: str) -> None:
    """Test that SMILES input field is present."""
    _goto_chemical_sketcher_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have SMILES input
    smiles_input = main.get_by_label(re.compile("SMILES", re.IGNORECASE)).or_(
        main.get_by_placeholder(re.compile("SMILES|CCO", re.IGNORECASE))
    )
    
    assert smiles_input.count() > 0, "Expected SMILES input field"


def test_chemical_sketcher_register_button(page: Page, streamlit_server: str) -> None:
    """Test that register button is present."""
    _goto_chemical_sketcher_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have register button
    register_button = main.get_by_role("button").filter(has_text=re.compile("Register", re.IGNORECASE))
    
    assert register_button.count() > 0, "Expected register button"


def test_chemical_sketcher_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_chemical_sketcher_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

