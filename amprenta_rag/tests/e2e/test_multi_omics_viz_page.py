"""E2E tests for Multi-Omics Integration dashboard page."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Multi-Omics%20Integration")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)


def test_multi_omics_integration_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Multi-Omics Integration page loads successfully."""
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Multi-Omics Integration"))).to_be_visible(timeout=20000)


def test_tabs_exist(page: Page, streamlit_server: str) -> None:
    """Test that all 6 tabs exist on the Multi-Omics Integration page."""
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Check all 6 tabs exist
    expect(main.locator('button[data-baseweb="tab"]', has_text="Setup + Run").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Variance").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Factor Scatter").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Loadings").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Alluvial").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="UpSet").first).to_be_visible(timeout=10000)


def test_alluvial_tab_exists(page: Page, streamlit_server: str) -> None:
    """Test that the Alluvial tab contains expected elements."""
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click on Alluvial tab
    alluvial_tab = main.locator('button[data-baseweb="tab"]', has_text="Alluvial").first
    alluvial_tab.click()
    
    # Wait for tab content to load
    page.wait_for_timeout(1000)
    
    # Check for key elements in the Alluvial tab
    expect(page.get_by_text("Alluvial (Dataset Flow)")).to_be_visible(timeout=5000)
    
    # Look for dataset selection multiselect (may show "No datasets available" if no data)
    # We check for either the multiselect or the "No datasets available" message
    dataset_selector_or_message = page.locator('text="Select datasets"').or_(page.locator('text="No datasets available"'))
    expect(dataset_selector_or_message.first).to_be_visible(timeout=5000)
    
    # Check for "Top N features" slider (should be visible regardless of data)
    expect(page.get_by_text("Top N features")).to_be_visible(timeout=5000)
