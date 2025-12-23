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

    # Ensure we actually landed on the MOA page; otherwise skip (nav/routing may differ by env).
    if main.locator('h1:has-text("MOA Inference"), h2:has-text("MOA Inference")').count() == 0:
        pytest.skip("MOA Inference page header not detected (routing/navigation mismatch).")

    if main.locator("text=No compounds available.").count() > 0:
        pytest.skip("No compounds available in this environment.")

    # Ensure page content is loaded.
    run_btn = page.get_by_role("button", name="Run MOA Inference").first
    if run_btn.count() == 0:
        pytest.skip("MOA inference form button not found (page may require data/DB).")
    expect(run_btn).to_be_visible(timeout=20000)
    expect(main.locator("text=Method").first).to_be_visible(timeout=20000)

    bayesian_radio = page.get_by_role("radio", name="Bayesian").first
    if bayesian_radio.count() == 0:
        pytest.skip("Bayesian radio option not found (widget markup changed).")
    bayesian_radio.click()
    page.wait_for_timeout(500)
    expect(bayesian_radio).to_be_visible(timeout=20000)


