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

    # Check for page heading - using get_by_text for partial match
    expect(main.get_by_text("Sphingolipid Pathway Imbalance").first).to_be_visible(timeout=20000)

    # Check for caption text that's always rendered (regardless of data state)
    # This proves the page fully loaded and rendered its content
    caption_or_datasets = main.get_by_text("Heuristic scoring").or_(
        main.get_by_text("Select datasets")
    ).or_(
        main.get_by_text("No datasets available.")
    )
    expect(caption_or_datasets.first).to_be_visible(timeout=15000)


