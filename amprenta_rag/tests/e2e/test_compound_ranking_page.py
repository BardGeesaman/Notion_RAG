"""E2E tests for Compound Ranking dashboard page."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=Compound%20Ranking")
    try:
        page.wait_for_load_state("domcontentloaded")
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)


def test_compound_ranking_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Compound Ranking page loads successfully."""
    _goto(page, streamlit_server)
    # Try multiple selectors to find the page content
    heading_selector = page.get_by_role("heading", name=re.compile(r"Compound Ranking"))
    text_selector = page.get_by_text("Compound Ranking")
    caption_selector = page.get_by_text("Multi-objective compound optimization and ranking")
    
    # Check if any of these elements are visible
    try:
        expect(heading_selector).to_be_visible(timeout=10000)
    except AssertionError:
        try:
            expect(text_selector.first).to_be_visible(timeout=10000)
        except AssertionError:
            expect(caption_selector).to_be_visible(timeout=10000)


# Note: Duplicate test functions removed - functionality covered by test_compound_ranking_page_loads
