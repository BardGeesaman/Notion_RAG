"""
E2E tests for optimistic locking conflict handling.

These tests exercise the Streamlit dashboard's concurrent-edit UX by opening the
same experiment in two browser contexts and verifying that:
- a stale save shows a conflict message,
- reload clears the conflict and allows saving,
- sequential saves in one session succeed.
"""

from __future__ import annotations

import time

import pytest
from playwright.sync_api import Browser, Page, expect

pytestmark = pytest.mark.requires_server


def _goto_experiment_edit_design(page: Page, base_url: str) -> None:
    """Navigate to Experiments -> Edit Design tab."""
    page.goto(base_url)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(5000)
    
    # Click Experiments button - may need multiple clicks for Streamlit
    experiments_btn = page.locator('button:has-text("Experiments")').first
    experiments_btn.wait_for(state="visible", timeout=10000)
    experiments_btn.click()
    page.wait_for_timeout(2000)
    
    # Check if page changed, if not click again
    if page.locator("text=Browse Experiments").count() == 0:
        experiments_btn.click()
        page.wait_for_timeout(3000)
    
    # Screenshot for debugging
    page.screenshot(path="/tmp/e2e_after_experiments.png")
    
    # Now look for Edit Design tab
    edit_tab = page.locator('button:has-text("Edit Design")').first
    if edit_tab.count() > 0:
        edit_tab.click()
        page.wait_for_timeout(2000)
    else:
        raise Exception("Edit Design tab not found - navigation may have failed")

    # Check if experiments are available
    no_experiments = page.locator("text=No experiments available to edit.")
    select_experiment = page.locator("text=Select experiment").first
    
    has_experiments = select_experiment.count() > 0
    has_no_data_msg = no_experiments.count() > 0
    
    # Assert that the page loaded properly
    assert has_experiments or has_no_data_msg, "Expected either experiment selector or 'No experiments available' message"
    
    # If no experiments, raise exception so test fails with clear message
    if has_no_data_msg:
        raise Exception("No experiments available to edit - test requires seeded data")
    
    expect(select_experiment).to_be_visible()


def _set_sample_groups_json(page: Page, payload: str) -> None:
    """Set the Sample groups JSON text area to a valid JSON string."""
    area = page.locator('textarea[aria-label="Sample groups (JSON)"]').first
    expect(area).to_be_visible()
    area.fill(payload)


def _click_save_design(page: Page) -> None:
    btn = page.get_by_role("button", name="Save design").first
    expect(btn).to_be_visible()
    btn.click()
    page.wait_for_timeout(1500)


class TestOptimisticLocking:
    def test_concurrent_experiment_edit_shows_conflict(self, page: Page, browser: Browser, streamlit_server: str):
        """Two users editing same experiment - second save shows conflict."""
        _goto_experiment_edit_design(page, streamlit_server)

        ctx2 = browser.new_context()
        page2 = ctx2.new_page()
        try:
            _goto_experiment_edit_design(page2, streamlit_server)

            t = int(time.time())
            _set_sample_groups_json(page, f'{{\"_e2e\": \"first-{t}\"}}')
            _click_save_design(page)
            expect(page.locator("text=Saved successfully!").first).to_be_visible()

            _set_sample_groups_json(page2, f'{{\"_e2e\": \"second-{t}\"}}')
            _click_save_design(page2)
            expect(page2.locator("text=⚠️ Conflict").first).to_be_visible()
        finally:
            ctx2.close()

    def test_reload_after_conflict_shows_updated_version(self, page: Page, browser: Browser, streamlit_server: str):
        """After conflict, reload button refreshes data and allows saving."""
        _goto_experiment_edit_design(page, streamlit_server)

        ctx2 = browser.new_context()
        page2 = ctx2.new_page()
        try:
            _goto_experiment_edit_design(page2, streamlit_server)

            t = int(time.time())
            first_payload = f'{{\"_e2e\": \"first-{t}\"}}'
            second_payload = f'{{\"_e2e\": \"second-{t}\"}}'

            _set_sample_groups_json(page, first_payload)
            _click_save_design(page)
            expect(page.locator("text=Saved successfully!").first).to_be_visible()

            _set_sample_groups_json(page2, second_payload)
            _click_save_design(page2)
            expect(page2.locator("text=⚠️ Conflict").first).to_be_visible()

            # Reload latest (clears conflict, refreshes internal version)
            reload_btn = page2.get_by_role("button", name="Reload Latest").first
            expect(reload_btn).to_be_visible()
            reload_btn.click()
            
            # Wait for page to fully reload after st.rerun()
            page2.wait_for_load_state("networkidle")
            page2.wait_for_timeout(2000)
            
            # Verify form is ready (sample groups textarea visible)
            ta = page2.locator('textarea[aria-label="Sample groups (JSON)"]')
            expect(ta).to_be_visible(timeout=10000)

            expect(page2.locator("text=⚠️ Conflict").first).to_have_count(0)

            # Now saving should succeed (no conflict) because the version was refreshed.
            _set_sample_groups_json(page2, f'{{\"_e2e\": \"after-reload-{t}\"}}')
            _click_save_design(page2)
            page2.screenshot(path="/tmp/e2e_after_reload_save.png")
            expect(page2.locator("text=Saved successfully!").first).to_be_visible(timeout=10000)
        finally:
            ctx2.close()

    def test_sequential_edits_succeed(self, page: Page, streamlit_server: str):
        """Sequential edits in one session succeed (no conflict)."""
        _goto_experiment_edit_design(page, streamlit_server)

        t = int(time.time())
        _set_sample_groups_json(page, f'{{\"_e2e\": \"seq-1-{t}\"}}')
        _click_save_design(page)
        expect(page.locator("text=Saved successfully!").first).to_be_visible()

        # No conflict should be shown in a single-session sequential save.
        expect(page.locator("text=⚠️ Conflict").first).to_have_count(0)

        _set_sample_groups_json(page, f'{{\"_e2e\": \"seq-2-{t}\"}}')
        _click_save_design(page)
        expect(page.locator("text=Saved successfully!").first).to_be_visible()
        expect(page.locator("text=⚠️ Conflict").first).to_have_count(0)


