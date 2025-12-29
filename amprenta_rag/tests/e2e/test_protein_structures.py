from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Protein%20Structures")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_protein_structures_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"Protein Structures"))).to_be_visible(timeout=20000)


def test_fetch_form_elements_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first

    expect(main.locator("text=Fetch structure").first).to_be_visible(timeout=10000)
    expect(main.locator("text=Source").first).to_be_visible(timeout=10000)
    expect(page.locator('input[aria-label="PDB ID"]').first).to_be_visible(timeout=10000)
    expect(page.locator('input[aria-label="UniProt ID"]').first).to_be_visible(timeout=10000)
    expect(page.get_by_role("button", name="Fetch").first).to_be_visible(timeout=10000)


def test_filter_dropdowns_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first

    expect(main.locator("text=Filter source").first).to_be_visible(timeout=10000)
    expect(main.locator("text=Filter prep status").first).to_be_visible(timeout=10000)


