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
    
    # Use label-based approach - check for visible "Use Graph Boost" text
    # This is more reliable across Streamlit versions than input[aria-label]
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for the checkbox label text (Streamlit renders checkbox with label)
    graph_boost_label = main.get_by_text("Use Graph Boost")
    expect(graph_boost_label.first).to_be_visible(timeout=20000)


