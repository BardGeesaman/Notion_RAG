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


def test_batch_correction_page_loads(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "Batch%20Correction")
    main = _main_container(page)
    expect(main.locator("text=Batch Effect Correction").first).to_be_visible(timeout=20000)


def test_batch_selection_and_assignment_ui(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "Batch%20Correction")
    main = _main_container(page)

    # If there are no datasets in DB, page displays a message and returns.
    if main.locator("text=No datasets available.").count() > 0:
        pytest.skip("No datasets available to exercise selection widgets.")

    expect(main.locator("text=Select datasets").first).to_be_visible(timeout=20000)

    # Try selecting the first dataset option in multiselect.
    ms = page.locator('input[aria-label="Select datasets"]')
    if ms.count() == 0:
        pytest.skip("Multiselect widget not found (Streamlit markup changed).")

    ms.first.click()
    page.wait_for_timeout(500)
    # pick first option if any
    opt = page.locator('[role="option"]').first
    if opt.count() == 0:
        pytest.skip("No multiselect options found to select.")
    opt.click()
    page.wait_for_timeout(500)

    expect(main.locator("text=Assign batches").first).to_be_visible(timeout=20000)


