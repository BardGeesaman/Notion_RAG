from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Compound-Target%20Network")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)


def test_compound_target_network_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Compound-Target Network"))).to_be_visible(timeout=20000)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.locator('button[data-baseweb="tab"]', has_text="Build Network").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Visualize").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Details").first).to_be_visible(timeout=10000)


