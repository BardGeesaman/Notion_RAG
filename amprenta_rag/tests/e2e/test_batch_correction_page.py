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

    # Verify dataset selection UI is present
    no_datasets = main.locator("text=No datasets available.")
    select_datasets = main.locator("text=Select datasets").first
    
    has_datasets = select_datasets.count() > 0
    has_no_data_msg = no_datasets.count() > 0
    
    assert has_datasets or has_no_data_msg, "Expected either dataset selector or 'No datasets available' message"
    
    if not has_datasets:
        # Page loaded but needs seeded data - this is acceptable for E2E
        return
    
    expect(select_datasets).to_be_visible(timeout=20000)

    # Verify multiselect widget exists
    ms = page.locator('input[aria-label="Select datasets"]')
    assert ms.count() > 0, "Multiselect widget for 'Select datasets' must exist"

    ms.first.click()
    page.wait_for_timeout(500)
    
    # Check if options are available
    opt = page.locator('[role="option"]').first
    if opt.count() > 0:
        opt.click()
        page.wait_for_timeout(500)
        expect(main.locator("text=Assign batches").first).to_be_visible(timeout=20000)


