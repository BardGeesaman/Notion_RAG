"""
E2E tests for the Nextflow Orchestrator dashboard page.

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


def _goto_nextflow_orchestrator(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Nextflow%20Orchestrator")
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)


class TestNextflowOrchestratorE2E:
    def test_nextflow_page_loads(self, page: Page, streamlit_server: str):
        """Test that Nextflow Orchestrator page loads with header."""
        _goto_nextflow_orchestrator(page, streamlit_server)

        main = _main_container(page)
        expect(
            main.locator('h1:has-text("Nextflow Orchestrator"), h2:has-text("Nextflow Orchestrator")').first
        ).to_be_visible(timeout=20000)

    def test_nextflow_has_tabs(self, page: Page, streamlit_server: str):
        """Test that all three tabs are present."""
        _goto_nextflow_orchestrator(page, streamlit_server)

        expect(page.get_by_role("tab", name="Run Pipeline")).to_be_visible(timeout=10000)
        expect(page.get_by_role("tab", name="Job Monitor")).to_be_visible(timeout=10000)
        expect(page.get_by_role("tab", name="Results")).to_be_visible(timeout=10000)


