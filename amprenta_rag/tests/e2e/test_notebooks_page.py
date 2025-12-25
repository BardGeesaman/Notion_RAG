from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Notebooks")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_notebooks_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.get_by_test_id("stMainBlockContainer").first
    expect(main.get_by_role("heading", name=re.compile(r"^Notebooks$"))).to_be_visible(timeout=20000)


def test_open_as_dashboard_button_exists(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    # At least one card should show "Open as Dashboard" if templates exist.
    expect(main.locator("text=Open as Dashboard").first).to_be_visible(timeout=20000)


def test_notebook_card_displays_metadata(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    # Cards should show basic metadata and action buttons.
    expect(main.locator("text=Open in Jupyter").first).to_be_visible(timeout=20000)
    expect(main.locator("text=View Source").first).to_be_visible(timeout=20000)


