from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Pathway%20Map%20Viewer")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)


def test_pathway_map_viewer_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Pathway Map Viewer"))).to_be_visible(timeout=20000)


def test_tabs_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.locator('button[data-baseweb="tab"]', has_text="Browse Pathways").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Pathway Map").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Data Overlay").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Multi-Pathway").first).to_be_visible(timeout=10000)


