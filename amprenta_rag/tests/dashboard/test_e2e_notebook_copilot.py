"""
E2E tests for the Notebook Co-Pilot dashboard page.

These tests validate that the page renders and that the copilot UI can generate
code and expose a copy button (when the environment supports generation).
"""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto_notebook_copilot(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Notebook%20Co-Pilot")
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)


def _skip_if_no_datasets(page: Page) -> None:
    # The selector displays a placeholder option when none exist.
    if page.locator("text=(no datasets found)").count() > 0:
        pytest.skip("No datasets available in this environment.")


class TestNotebookCopilotE2E:
    def test_copilot_page_loads(self, page: Page, streamlit_server: str):
        """Verify Notebook Co-Pilot page renders with key elements."""
        _goto_notebook_copilot(page, streamlit_server)

        expect(page.locator("text=Notebook Co-Pilot").first).to_be_visible(timeout=10000)
        expect(page.get_by_role("button", name="Load Dataset").first).to_be_visible()
        expect(page.get_by_role("button", name="HTS QC").first).to_be_visible()
        expect(page.get_by_role("button", name="Dose-response").first).to_be_visible()
        expect(page.get_by_role("button", name="Publish").first).to_be_visible()
        expect(page.get_by_role("button", name="Generate code").first).to_be_visible()

    def test_load_dataset_action(self, page: Page, streamlit_server: str):
        """Click Load Dataset and verify code preview can appear after generation."""
        _goto_notebook_copilot(page, streamlit_server)

        page.get_by_role("button", name="Load Dataset").first.click()
        page.wait_for_timeout(500)

        # Ensure the dataset selector is present; if none exist, skip.
        expect(page.locator("text=Dataset").first).to_be_visible(timeout=10000)
        _skip_if_no_datasets(page)

        page.get_by_role("button", name="Generate code").first.click()
        page.wait_for_timeout(3000)

        if page.locator("text=Copilot failed:").count() > 0:
            pytest.skip("Copilot generation unavailable (missing LLM credentials/config).")

        expect(page.locator("text=Generated code").first).to_be_visible(timeout=20000)
        expect(page.locator("css=pre").first).to_be_visible(timeout=20000)

    def test_copy_button_visible(self, page: Page, streamlit_server: str):
        """Verify copy button shows after successful code generation."""
        _goto_notebook_copilot(page, streamlit_server)

        page.get_by_role("button", name="Load Dataset").first.click()
        page.wait_for_timeout(500)
        _skip_if_no_datasets(page)

        page.get_by_role("button", name="Generate code").first.click()
        page.wait_for_timeout(3000)

        if page.locator("text=Copilot failed:").count() > 0:
            pytest.skip("Copilot generation unavailable (missing LLM credentials/config).")

        expect(page.get_by_role("button", name="Copy").first).to_be_visible(timeout=20000)


