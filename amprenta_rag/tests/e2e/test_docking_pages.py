from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str, page_name: str) -> None:
    page.goto(f"{base_url}/?page={page_name}")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_docking_runs_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server, "Docking%20Runs")
    expect(page.get_by_role("heading", name=re.compile(r"Docking Runs"))).to_be_visible(timeout=20000)


def test_docking_triage_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server, "Docking%20Triage")
    expect(page.get_by_role("heading", name=re.compile(r"Docking Triage"))).to_be_visible(timeout=20000)


