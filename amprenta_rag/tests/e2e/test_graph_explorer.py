from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Graph%20Explorer")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        # Streamlit can keep long-lived connections; fall back to a shorter readiness heuristic.
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_graph_explorer_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Graph Explorer"))).to_be_visible(timeout=20000)


def test_traverse_form_elements_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first

    # Traverse section elements
    expect(main.locator("text=Entity type").first).to_be_visible(timeout=10000)
    expect(page.locator('input[aria-label="Entity UUID"]').first).to_be_visible(timeout=10000)
    expect(main.locator("text=k-hop").first).to_be_visible(timeout=10000)
    expect(page.get_by_role("button", name="Traverse").first).to_be_visible(timeout=10000)


def test_path_form_elements_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first

    expect(main.locator("text=Shortest path").first).to_be_visible(timeout=10000)
    expect(main.locator("text=Source type").first).to_be_visible(timeout=10000)
    expect(page.locator('input[aria-label="Source UUID"]').first).to_be_visible(timeout=10000)
    expect(main.locator("text=Target type").first).to_be_visible(timeout=10000)
    expect(page.locator('input[aria-label="Target UUID"]').first).to_be_visible(timeout=10000)
    expect(page.get_by_role("button", name=re.compile(r"Find path")).first).to_be_visible(timeout=10000)


def test_analytics_checkboxes_exist(page: Page, streamlit_server: str) -> None:
    """Verify analytics toggles exist on Graph Explorer page."""
    _goto(page, streamlit_server)
    expect(page.get_by_text("Size by Degree").first).to_be_visible(timeout=10000)
    expect(page.get_by_text("Color by Community").first).to_be_visible(timeout=10000)


