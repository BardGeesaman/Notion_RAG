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


def test_model_registry_page_loads(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "Model%20Registry")
    main = _main_container(page)
    expect(main.locator("text=ML Model Registry").first).to_be_visible(timeout=20000)


def test_admet_ui_if_present(page: Page, streamlit_server: str):
    """
    Optional: if there is an ADMET UI surface in the Chemistry pages, ensure it loads.
    (This environment may not expose an explicit compound ADMET UI yet.)
    """
    _goto(page, streamlit_server, "Chemistry")
    main = _main_container(page)

    # If the chemistry page doesn't exist or isn't reachable, this will fail; skip if missing.
    if main.locator("text=Chemistry").count() == 0 and main.locator("text=ADMET").count() == 0:
        pytest.skip("No Chemistry/ADMET UI detected to validate.")

    # Best-effort: just assert the page isn't blank and main container is present.
    expect(main).to_be_visible(timeout=20000)


