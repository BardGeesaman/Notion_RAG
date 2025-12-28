from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _main_container(page: Page):
    main = page.locator('[data-testid="stMainBlockContainer"]')
    if main.count() > 0:
        return main.first
    return page.locator('[data-testid="stAppViewContainer"]').first


def _goto(page: Page, base_url: str, page_name: str) -> None:
    page.goto(f"{base_url}/?page={page_name}")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_sphingolipid_imbalance_page_loads(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "Sphingolipid%20Imbalance")
    main = _main_container(page)

    expect(main.locator("text=Sphingolipid Pathway Imbalance").first).to_be_visible(timeout=20000)

    # Verify dataset selection UI is present
    # If no datasets message appears, the page loaded but needs seeded data
    no_datasets = main.locator("text=No datasets available.")
    select_datasets = main.locator("text=Select datasets").first
    
    # Assert either the selection widget OR a clear message about needing data
    has_datasets = select_datasets.count() > 0
    has_no_data_msg = no_datasets.count() > 0
    
    assert has_datasets or has_no_data_msg, "Expected either dataset selector or 'No datasets available' message"
    
    # If datasets are available, verify the selector is visible
    if has_datasets:
        expect(select_datasets).to_be_visible(timeout=20000)


