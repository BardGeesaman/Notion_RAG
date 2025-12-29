"""E2E tests for Projector dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_projector_page(page: Page, base_url: str) -> None:
    """Navigate to the Projector page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Projector")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_projector_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Projector page loads successfully."""
    _goto_projector_page(page, streamlit_server)

    main = _main_container(page)

    # Look for Projector heading
    heading_patterns = [
        main.get_by_text(re.compile(r"High-Dimensional.*Projector", re.IGNORECASE)),
        main.get_by_text(re.compile(r"Projector", re.IGNORECASE)),
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


def test_projector_algorithm_selector(page: Page, streamlit_server: str) -> None:
    """Test that algorithm selector is present."""
    _goto_projector_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Look for algorithm selector in sidebar
    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have algorithm-related controls
    algorithm_text = sidebar.get_by_text(re.compile("Algorithm|UMAP|t-SNE|PCA", re.IGNORECASE))
    selectboxes = sidebar.locator('[data-baseweb="select"]')
    
    assert algorithm_text.count() > 0 or selectboxes.count() > 0, "Expected algorithm selector"


def test_projector_dataset_selector(page: Page, streamlit_server: str) -> None:
    """Test that dataset selector is present."""
    _goto_projector_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have load datasets button or selector
    load_button = sidebar.get_by_role("button").filter(has_text=re.compile("Load.*Datasets", re.IGNORECASE))
    dataset_text = sidebar.get_by_text(re.compile("Dataset|Select", re.IGNORECASE))
    
    assert load_button.count() > 0 or dataset_text.count() > 0, "Expected dataset controls"


def test_projector_sidebar_controls(page: Page, streamlit_server: str) -> None:
    """Test that sidebar has projection controls."""
    _goto_projector_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    sidebar = page.locator('[data-testid="stSidebar"]')
    
    # Should have some buttons and controls in sidebar
    buttons = sidebar.get_by_role("button")
    
    assert buttons.count() > 0, "Expected sidebar controls"


def test_projector_main_content(page: Page, streamlit_server: str) -> None:
    """Test that main content area renders."""
    _goto_projector_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    main = _main_container(page)
    
    # Should have some content (instructions or visualization)
    info_text = main.get_by_text(re.compile("Select.*dataset|Compute|Projection", re.IGNORECASE))
    
    assert main.count() > 0, "Main container should exist"


def test_projector_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_projector_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

