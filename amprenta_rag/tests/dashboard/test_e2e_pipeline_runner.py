"""
E2E tests for the Pipeline Runner dashboard page.

Tests validate page loads, tabs present, and basic UI structure.
"""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _main_container(page: Page):
    """Return locator for Streamlit main content region."""
    main = page.locator('[data-testid="stMainBlockContainer"]')
    if main.count() > 0:
        return main.first
    return page.locator('[data-testid="stAppViewContainer"]').first


def _goto_pipeline_runner(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Pipeline%20Runner")
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)


class TestPipelineRunnerE2E:
    def test_pipeline_runner_page_loads(self, page: Page, streamlit_server: str):
        """Test that Pipeline Runner page loads with header."""
        _goto_pipeline_runner(page, streamlit_server)

        main = _main_container(page)
        expect(
            main.locator('h1:has-text("Pipeline Runner"), h2:has-text("Pipeline Runner")').first
        ).to_be_visible(timeout=20000)

    def test_pipeline_runner_has_tabs(self, page: Page, streamlit_server: str):
        """Test that all three tabs are present."""
        _goto_pipeline_runner(page, streamlit_server)

        expect(page.get_by_role("tab", name="Run Pipeline")).to_be_visible(timeout=10000)
        expect(page.get_by_role("tab", name="Manage Indices")).to_be_visible(timeout=10000)
        expect(page.get_by_role("tab", name="Job History")).to_be_visible(timeout=10000)

    def test_run_pipeline_tab_structure(self, page: Page, streamlit_server: str):
        """Test Run Pipeline tab has expected elements."""
        _goto_pipeline_runner(page, streamlit_server)

        # Click Run Pipeline tab (should be default)
        page.get_by_role("tab", name="Run Pipeline").click()
        page.wait_for_timeout(1000)

        main = _main_container(page)
        # Should have input section
        expect(main.locator('text=Select Input').first).to_be_visible(timeout=10000)
        # Should have tool selector
        expect(main.locator('label:has-text("Tool")').first).to_be_visible(timeout=10000)

    def test_manage_indices_tab_structure(self, page: Page, streamlit_server: str):
        """Test Manage Indices tab has expected elements."""
        _goto_pipeline_runner(page, streamlit_server)

        page.get_by_role("tab", name="Manage Indices").click()
        page.wait_for_timeout(1000)

        main = _main_container(page)
        # Should have upload/register form action
        expect(main.get_by_role("button", name="Register Index").first).to_be_visible(timeout=10000)
        # Should have organism field
        expect(main.locator('label:has-text("Organism")').first).to_be_visible(timeout=10000)


