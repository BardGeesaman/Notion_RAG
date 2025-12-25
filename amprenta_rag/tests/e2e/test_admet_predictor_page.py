from __future__ import annotations

import re

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


def _goto(page: Page, base_url: str) -> None:
    page.goto(f"{base_url}/?page=ADMET%20Predictor")
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2500)


def test_admet_predictor_page_loads(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    expect(page.get_by_role("heading", name=re.compile(r"ADMET Predictor"))).to_be_visible(timeout=20000)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    # Tabs are rendered via BaseWeb.
    expect(main.locator('button[data-baseweb="tab"]', has_text="Predict").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Calibration").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Model Info").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Explain").first).to_be_visible(timeout=10000)
    expect(main.locator('button[data-baseweb="tab"]', has_text="Global Importance").first).to_be_visible(timeout=10000)


def test_predict_tab_elements_exist(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_text("SMILES (one per line, max 100)").first).to_be_visible(timeout=10000)
    expect(main.locator('textarea[aria-label="SMILES (one per line, max 100)"]').first).to_be_visible(timeout=10000)
    expect(main.get_by_text("Endpoints").first).to_be_visible(timeout=10000)
    expect(main.get_by_text("Show uncertainty").first).to_be_visible(timeout=10000)
    expect(main.locator('button:not([data-baseweb="tab"])', has_text=re.compile(r"^Predict$")).first).to_be_visible(
        timeout=10000
    )


def test_predict_flow_with_valid_smiles(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    main.locator('textarea[aria-label="SMILES (one per line, max 100)"]').first.fill("CCO")
    page.keyboard.press("Tab")
    main.locator('button:not([data-baseweb="tab"])', has_text=re.compile(r"^Predict$")).first.click()
    page.wait_for_timeout(2500)

    # Either results appear (dataframe) OR we show a no-model/rdkit message; both still render Download CSV.
    expect(main.locator('[data-testid="stDownloadButton"]').first).to_be_visible(timeout=30000)


def test_predict_invalid_smiles_shows_error(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    main.locator('textarea[aria-label="SMILES (one per line, max 100)"]').first.fill("INVALID_SMILES")
    page.keyboard.press("Tab")
    main.locator('button:not([data-baseweb="tab"])', has_text=re.compile(r"^Predict$")).first.click()
    page.wait_for_timeout(2000)

    # Error text should be visible (either invalid smiles or missing RDKit in some environments).
    expect(
        main.locator('[data-testid="stAlert"]', has_text=re.compile(r"(Invalid SMILES|RDKit required)")).first
    ).to_be_visible(timeout=20000)


def test_model_info_tab_shows_registry(page: Page, streamlit_server: str) -> None:
    _goto(page, streamlit_server)
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    main.get_by_text("Model Info").first.click()
    page.wait_for_timeout(1000)
    # Either registry table appears or the fallback message.
    expect(
        main.get_by_text(re.compile(r"(Registered ADMET models|No models found|No admet_\\* models)")).first
    ).to_be_visible(timeout=20000)


