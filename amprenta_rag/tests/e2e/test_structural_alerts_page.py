from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Structural%20Alerts")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)


def test_structural_alerts_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Structural Alerts"))).to_be_visible(timeout=20000)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.locator('button[data-baseweb="tab"]', has_text="Single Compound").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Batch Check").first).to_be_visible(timeout=10000)


def test_single_compound_check(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    main.locator('input[aria-label="SMILES"]').first.fill("CCO")
    page.keyboard.press("Tab")
    main.locator('button:not([data-baseweb="tab"])', has_text=re.compile(r"^Check Alerts$")).first.click()
    page.wait_for_timeout(2000)
    expect(main.get_by_text("GREEN").first).to_be_visible(timeout=20000)



