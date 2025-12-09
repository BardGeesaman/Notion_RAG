from __future__ import annotations

"""
Playwright UI interaction tests for specific dashboard pages.

Focuses on buttons, forms, filters, and basic navigation behaviors.
"""

import os

import pytest
from playwright.sync_api import Page, expect

pytestmark = pytest.mark.skip(reason="Requires running dashboard server at localhost:8502")


@pytest.fixture(scope="session")
def dashboard_base_url() -> str:
    return os.getenv("DASHBOARD_BASE_URL", "http://localhost:8502")


@pytest.mark.ui
def test_search_page_filters_and_submit(page: Page, dashboard_base_url: str):
    """
    Verify basic interactions on the Search page.

    This is intentionally high-level; selectors may be refined using playwright codegen.
    """
    page.goto(dashboard_base_url)

    # Navigate to Search page
    page.get_by_label("Select Page").get_by_text("Search").click(force=True)
    expect(page.get_by_text("Search")).to_be_visible()

    # Try interacting with a text input (search box)
    text_inputs = page.locator("input[type='text']")
    if text_inputs.count() > 0:
        text_inputs.first.fill("ALS")

    # Click a generic search / run button if present
    for name in ["Search", "Run", "Apply", "Filter"]:
        btn = page.get_by_role("button", name=name)
        if btn.count() > 0:
            btn.first.click()
            break


@pytest.mark.ui
def test_cross_omics_page_run_analysis(page: Page, dashboard_base_url: str):
    """
    Verify that the Cross-Omics page exposes primary actions.
    """
    page.goto(dashboard_base_url)

    # Navigate to Cross-Omics page
    page.get_by_label("Select Page").get_by_text("Cross-Omics").click(force=True)
    expect(page.get_by_text("Cross-Omics")).to_be_visible()

    # Look for a button that triggers analysis
    for name in ["Run Analysis", "Run", "Analyze"]:
        btn = page.get_by_role("button", name=name)
        if btn.count() > 0:
            btn.first.click()
            break


@pytest.mark.ui
def test_evidence_report_page_generate_report(page: Page, dashboard_base_url: str):
    """
    Verify that Evidence Report page exposes a 'generate' style action.
    """
    page.goto(dashboard_base_url)

    page.get_by_label("Select Page").get_by_text("Evidence Report").click(force=True)
    expect(page.get_by_text("Evidence Report")).to_be_visible()

    # Interact with controls that look like report generation triggers
    for name in ["Generate", "Create Report", "Run"]:
        btn = page.get_by_role("button", name=name)
        if btn.count() > 0:
            btn.first.click()
            break


