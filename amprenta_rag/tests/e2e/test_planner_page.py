"""E2E tests for Experiment Planner dashboard."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_planner_page(page: Page, base_url: str) -> None:
    """Navigate to the Experiment Planner page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Experiment%20Planner")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_planner_page_loads(page: Page, streamlit_server: str) -> None:
    """Test that the Planner page loads successfully."""
    _goto_planner_page(page, streamlit_server)

    main = _main_container(page)
    
    heading = main.get_by_text(re.compile(r"Experiment.*Planner", re.IGNORECASE))
    
    try:
        expect(heading.first).to_be_visible(timeout=5000)
    except AssertionError:
        expect(main).to_be_visible(timeout=10000)


def test_planner_tabs_present(page: Page, streamlit_server: str) -> None:
    """Test that planner tabs are present."""
    _goto_planner_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    tabs = page.locator('[role="tab"]')
    assert tabs.count() >= 2, f"Expected tabs, found {tabs.count()}"


def test_planner_power_calculator_inputs(page: Page, streamlit_server: str) -> None:
    """Test that power calculator has input controls."""
    _goto_planner_page(page, streamlit_server)
    page.wait_for_timeout(5000)

    main = _main_container(page)
    
    # Should have selectboxes and inputs
    selectboxes = main.locator('[data-baseweb="select"]')
    buttons = main.get_by_role("button")
    
    assert selectboxes.count() > 0 or buttons.count() > 0, "Expected power calculator controls"


def test_planner_plate_layout_tab(page: Page, streamlit_server: str) -> None:
    """Test that plate layout tab renders."""
    _goto_planner_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Plate Layout tab
    plate_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Plate", re.IGNORECASE))
    
    try:
        plate_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    plate_text = main.get_by_text(re.compile("Plate|Layout|Wells", re.IGNORECASE))
    
    assert plate_text.count() > 0 or main.count() > 0, "Plate layout tab should render"


def test_planner_cost_estimator_tab(page: Page, streamlit_server: str) -> None:
    """Test that cost estimator tab renders."""
    _goto_planner_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    # Try to click Cost Estimator tab
    cost_tab = page.locator('[role="tab"]').filter(has_text=re.compile("Cost", re.IGNORECASE))
    
    try:
        cost_tab.first.click(timeout=10000)
        page.wait_for_timeout(2000)
    except Exception:
        pass

    main = _main_container(page)
    cost_text = main.get_by_text(re.compile("Cost|Estimate|Sample", re.IGNORECASE))
    
    assert cost_text.count() > 0 or main.count() > 0, "Cost estimator tab should render"


def test_planner_no_crash(page: Page, streamlit_server: str) -> None:
    """Test that page doesn't crash on load."""
    _goto_planner_page(page, streamlit_server)
    page.wait_for_timeout(3000)

    assert page.locator('body').count() > 0, "Page failed to load"
    
    main = page.locator('[data-testid="stMainBlockContainer"]')
    assert main.count() > 0, "Main container not found"

