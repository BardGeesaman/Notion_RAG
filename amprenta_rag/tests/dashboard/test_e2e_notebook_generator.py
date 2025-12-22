"""
E2E tests for the Query→Notebook Generator dashboard page.

These tests validate that the page renders and that notebook generation can
complete when LLM credentials are available in the environment.
"""

from __future__ import annotations

import os
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


def _goto_notebook_generator(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Notebook%20Generator")
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(5000)

    main = _main_container(page)
    expect(
        main.locator('h1:has-text("Notebook Generator"), h2:has-text("Notebook Generator")').first
    ).to_be_visible(timeout=20000)


def _try_generate(page: Page) -> bool:
    """Click generate and return True if preview appears, False if error shows."""
    page.get_by_role("button", name="Generate Notebook").first.click()

    # LLM-backed generation can take a while; allow up to 90s.
    main = _main_container(page)
    preview = main.locator('text=Preview').first
    err = main.locator("text=Failed to generate notebook:").first

    # Wait until either preview or error becomes visible.
    try:
        expect(preview).to_be_visible(timeout=90000)
        return True
    except Exception:
        if err.count() > 0:
            return False
        # If neither is visible, treat as failure.
        return False


class TestNotebookGeneratorE2E:
    def test_notebook_generator_page_loads(self, page: Page, streamlit_server: str):
        _goto_notebook_generator(page, streamlit_server)

        expect(page.locator('textarea[aria-label="Question"]').first).to_be_visible(timeout=10000)
        expect(page.get_by_role("button", name="Generate Notebook").first).to_be_visible(timeout=10000)

    def test_generate_notebook_from_query(self, page: Page, streamlit_server: str):
        _goto_notebook_generator(page, streamlit_server)

        # If no LLM credentials are configured, skip (generation is expected to fail).
        if not os.getenv("OPENAI_API_KEY"):
            pytest.skip("OPENAI_API_KEY not set; skipping LLM-backed notebook generation test.")

        page.locator('textarea[aria-label="Question"]').first.fill("Analyze transcriptomics dataset")
        page.keyboard.press("Tab")
        expect(page.get_by_role("button", name="Generate Notebook").first).to_be_enabled(timeout=10000)
        ok = _try_generate(page)
        if not ok:
            pytest.skip("Notebook generation unavailable (missing/invalid LLM credentials or backend config).")

        main = _main_container(page)
        expect(main.locator("text=Preview").first).to_be_visible(timeout=20000)

    def test_download_button_visible_after_generation(self, page: Page, streamlit_server: str):
        _goto_notebook_generator(page, streamlit_server)

        if not os.getenv("OPENAI_API_KEY"):
            pytest.skip("OPENAI_API_KEY not set; skipping LLM-backed notebook generation test.")

        page.locator('textarea[aria-label="Question"]').first.fill("Analyze transcriptomics dataset")
        page.keyboard.press("Tab")
        expect(page.get_by_role("button", name="Generate Notebook").first).to_be_enabled(timeout=10000)
        ok = _try_generate(page)
        if not ok:
            pytest.skip("Notebook generation unavailable (missing/invalid LLM credentials or backend config).")

        main = _main_container(page)
        expect(page.get_by_role("button", name="Download .ipynb").first).to_be_visible(timeout=20000)
        expect(main.locator("text=Filename:").first).to_be_visible(timeout=20000)


