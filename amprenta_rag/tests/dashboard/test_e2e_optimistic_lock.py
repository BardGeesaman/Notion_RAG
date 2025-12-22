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
    page.wait_for_timeout(3000)

    # Click Experiments button
    page.locator('button:has-text("Experiments")').first.click()
    
    # Wait for Experiments page to load (unique header)
    expect(page.locator("text=🔬 Experiments").first).to_be_visible(timeout=10000)
    page.wait_for_timeout(1000)

    # Click Edit Design tab
    page.locator('button:has-text("Edit Design")').first.click()
    
    # Wait for tab content
    page.wait_for_timeout(2000)

    # Skip if no experiments
    if page.locator("text=No experiments available to edit.").count() > 0:
        pytest.skip("No experiments available to edit in this environment.")

    expect(page.locator("text=Select experiment").first).to_be_visible()


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
            page2.wait_for_timeout(2500)

            expect(page2.locator("text=⚠️ Conflict").first).to_have_count(0)

            # Now saving should succeed (no conflict) because the version was refreshed.
            _set_sample_groups_json(page2, f'{{\"_e2e\": \"after-reload-{t}\"}}')
            _click_save_design(page2)
            expect(page2.locator("text=Saved successfully!").first).to_be_visible()
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


