from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=RAG%20Query")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_use_graph_boost_checkbox_exists(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    # Streamlit checkbox is typically rendered as an input with aria-label equal to the label.
    cb = page.locator('input[aria-label="Use Graph Boost"]').first
    if cb.count() == 0:
        # Fallback: locate label text
        cb = page.locator("text=Use Graph Boost").first
    expect(cb).to_be_visible(timeout=20000)


