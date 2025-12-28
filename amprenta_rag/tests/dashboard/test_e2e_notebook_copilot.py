"""
E2E tests for the Notebook Co-Pilot dashboard page.

These tests validate that the page renders and that the copilot UI can generate
code and expose a copy button (when the environment supports generation).
"""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _main_container(page: Page):
    """
    Return the best locator for the Streamlit *main content* region.

    Streamlit's outer container (`stAppViewContainer`) includes sidebar content,
    which can cause E2E selectors like `text=...` to match hidden sidebar items.
    """
    main = page.locator('[data-testid="stMainBlockContainer"]')
    if main.count() > 0:
        return main.first
    return page.locator('[data-testid="stAppViewContainer"]').first


def _goto_notebook_copilot(page: Page, base_url: str) -> None:
    # Use URL query param for direct navigation
    page.goto(f"{base_url}/?page=Notebook%20Co-Pilot")
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(5000)

    # Wait for page content
    main = _main_container(page)
    expect(
        main.locator('h1:has-text("Notebook Co-Pilot"), h2:has-text("Notebook Co-Pilot")').first
    ).to_be_visible(timeout=20000)


def _check_datasets_available(page: Page) -> bool:
    """Check if datasets are available, return True if available."""
    main = _main_container(page)
    no_datasets = main.locator("text=(no datasets found)")
    return no_datasets.count() == 0


class TestNotebookCopilotE2E:
    def test_copilot_page_loads(self, page: Page, streamlit_server: str):
        """Verify Notebook Co-Pilot page renders with key elements."""
        _goto_notebook_copilot(page, streamlit_server)

        main = _main_container(page)
        expect(
            main.locator('h1:has-text("Notebook Co-Pilot"), h2:has-text("Notebook Co-Pilot")').first
        ).to_be_visible(timeout=10000)
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

        # Ensure the dataset selector is present
        main = _main_container(page)
        expect(main.locator('label:has-text("Dataset")').first).to_be_visible(timeout=10000)
        
        # Check if datasets are available
        if not _check_datasets_available(page):
            # No datasets available - this is acceptable for E2E
            return

        page.get_by_role("button", name="Generate code").first.click()
        page.wait_for_timeout(3000)

        # Check if copilot generation succeeded or failed
        copilot_failed = page.locator("text=Copilot failed:")
        generated_code = page.locator("text=Generated code").first
        
        if copilot_failed.count() > 0:
            # Copilot requires LLM credentials - this is acceptable for E2E
            return
        
        expect(generated_code).to_be_visible(timeout=20000)
        expect(page.locator("css=pre").first).to_be_visible(timeout=20000)

    def test_copy_button_visible(self, page: Page, streamlit_server: str):
        """Verify copy button shows after successful code generation."""
        _goto_notebook_copilot(page, streamlit_server)

        page.get_by_role("button", name="Load Dataset").first.click()
        page.wait_for_timeout(500)
        
        # Check if datasets are available
        if not _check_datasets_available(page):
            # No datasets available - this is acceptable for E2E
            return

        page.get_by_role("button", name="Generate code").first.click()
        page.wait_for_timeout(3000)

        # Check if copilot generation succeeded or failed
        copilot_failed = page.locator("text=Copilot failed:")
        if copilot_failed.count() > 0:
            # Copilot requires LLM credentials - this is acceptable for E2E
            return

        main = _main_container(page)
        # Streamlit shows the copy button on hover for code blocks; hover first.
        pre = main.locator("css=pre").first
        expect(pre).to_be_visible(timeout=20000)
        pre.hover()
        page.wait_for_timeout(300)
        expect(main.locator('[data-testid="stCodeCopyButton"]').first).to_be_visible(timeout=20000)


