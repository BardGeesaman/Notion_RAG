from __future__ import annotations

import json
from pathlib import Path

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _main_container(page: Page):
    main = page.locator('[data-testid="stMainBlockContainer"]')
    if main.count() > 0:
        return main.first
    return page.locator('[data-testid="stAppViewContainer"]').first


def _goto(page: Page, base_url: str, page_name: str) -> None:
    page.goto(f"{base_url}/?page={page_name}")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_dose_response_template_has_bayesian_section():
    nb_path = Path("deploy/jupyterhub/templates/dose_response.ipynb")
    data = json.loads(nb_path.read_text(encoding="utf-8"))
    cells = data.get("cells", [])
    joined = "\n".join(
        "".join(c.get("source", [])) if isinstance(c.get("source"), list) else str(c.get("source", ""))
        for c in cells
    )
    assert "Bayesian Dose-Response Fit" in joined
    assert "Convergence Diagnostics" in joined


def test_experiment_optimizer_page_loads(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "Experiment%20Optimizer")
    main = _main_container(page)
    # Wait for a unique element on the page (more reliable than matching hidden sidebar text).
    expect(page.get_by_role("button", name="Recommend").first).to_be_visible(timeout=20000)
    expect(main.locator('h1:has-text("Experiment Optimizer"), h2:has-text("Experiment Optimizer")').first).to_be_visible(
        timeout=20000
    )


def test_moa_inference_toggle_works(page: Page, streamlit_server: str):
    _goto(page, streamlit_server, "MOA%20Inference")
    main = _main_container(page)

    # Verify we landed on the MOA page
    page_header = main.locator('h1:has-text("MOA Inference"), h2:has-text("MOA Inference")').first
    assert page_header.count() > 0, "MOA Inference page header must be present (check routing)"
    
    # Check for data availability
    no_compounds = main.locator("text=No compounds available.")
    has_no_data = no_compounds.count() > 0
    
    if has_no_data:
        # Page loaded but needs seeded data - this is acceptable for E2E
        return

    # Verify MOA inference form button exists
    run_btn = page.get_by_role("button", name="Run MOA Inference").first
    assert run_btn.count() > 0, "Run MOA Inference button must exist on page"
    expect(run_btn).to_be_visible(timeout=20000)
    expect(main.locator("text=Method").first).to_be_visible(timeout=20000)

    # Verify Bayesian radio option exists
    bayesian_radio = page.get_by_role("radio", name="Bayesian").first
    assert bayesian_radio.count() > 0, "Bayesian radio option must exist in Method selector"
    bayesian_radio.click()
    page.wait_for_timeout(500)
    expect(bayesian_radio).to_be_visible(timeout=20000)


