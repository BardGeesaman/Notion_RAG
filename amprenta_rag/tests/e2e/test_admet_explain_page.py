from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytest.importorskip("playwright")

pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=ADMET%20Predictor")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_explain_tab_exists(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.locator('button[data-baseweb="tab"]', has_text="Explain").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Global Importance").first).to_be_visible(timeout=10000)


def test_explain_flow_shows_waterfall(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    main.locator('button[data-baseweb="tab"]', has_text="Explain").first.click()
    page.wait_for_timeout(750)

    # Fill SMILES + click Explain
    main.locator('input[aria-label="SMILES (single compound)"]').first.fill("CCO")
    page.keyboard.press("Tab")
    main.locator('button:not([data-baseweb="tab"])', has_text=re.compile(r"^Explain$")).first.click()
    page.wait_for_timeout(2000)

    # Waterfall is rendered via Plotly; selectors vary by Streamlit version, so use a robust class-based check.
    expect(main.locator(".js-plotly-plot").first).to_be_visible(timeout=20000)
    expect(main.get_by_text(re.compile(r"Download SHAP JSON")).first).to_be_visible(timeout=20000)


