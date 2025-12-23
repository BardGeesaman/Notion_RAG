from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=CRISPR%20Analysis")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_crispr_analysis_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"CRISPR Analysis"))).to_be_visible(timeout=20000)


