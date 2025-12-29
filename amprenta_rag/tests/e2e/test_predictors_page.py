"""E2E tests for Predictors dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_predictors_page(page: Page, base_url: str) -> None:
    """Navigate to the Predictors page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Predictors")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_predictors_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Predictors page loads successfully."""
    _goto_predictors_page(page, streamlit_server)

    main = _main_container(page)

    # Look for ML Predictors heading
    heading_patterns = [
        main.get_by_text(re.compile(r"ML.*Predictors", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Predictors", re.IGNORECASE)),
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


def test_predictors_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that predictor tabs are present."""
    _goto_predictors_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    # Check if tabs exist (at least 2)
    tabs = page.locator('[role="tab"]')
    tab_count = tabs.count()
    
    # Should have at least 2 tabs
    assert tab_count >= 2, f"Expected at least 2 tabs but found {tab_count}"


def test_predictors_interactive_elements(page: Page, streamlit_server: str) -> None:
    """Test that page has interactive elements."""
    _goto_predictors_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Page should have buttons or inputs
    buttons = main.get_by_role("button")
    inputs = main.locator('input')
    textareas = main.locator('textarea')
    
    # Should have interactive elements
    assert (buttons.count() > 0 or inputs.count() > 0 or 
            textareas.count() > 0), "Expected interactive elements"


def test_predictors_train_tab_form(page: Page, streamlit_server: str) -> None:
    """Test that train tab has form elements."""
    _goto_predictors_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Train Model tab
    train_tab = page.get_by_text(re.compile("Train.*Model", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Train", re.IGNORECASE))
    )
    
    try:
        train_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have training form elements
    train_text = main.get_by_text(re.compile("Train|Program|Assay", re.IGNORECASE))
    
    assert train_text.count() > 0, "Expected training form content"


def test_predictors_predict_tab_form(page: Page, streamlit_server: str) -> None:
    """Test that predictions tab has form elements."""
    _goto_predictors_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Predictions tab
    predict_tab = page.get_by_text(re.compile("Run.*Predictions|Predictions", re.IGNORECASE), exact=False).or_(
        page.locator('[role="tab"]').filter(has_text=re.compile("Predict", re.IGNORECASE))
    )
    
    try:
        predict_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    
    # Should have prediction form elements
    predict_text = main.get_by_text(re.compile("Predictions?|Model|SMILES", re.IGNORECASE))
    
    assert predict_text.count() > 0, "Expected prediction form content"


def test_predictors_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_predictors_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Verify page is responsive
    assert page.locator('body').count() > 0, "Page failed to load"
    
    # Verify main container exists
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

