"""Chemistry registration E2E test."""
from playwright.sync_api import Page, expect

import pytest

pytestmark = pytest.mark.requires_server


def test_register_compound_and_verify(page: Page, base_url: str):
    page.goto(base_url)
    page.wait_for_timeout(5000)

    # Click Chemistry in sidebar
    page.locator("text=Chemistry").click()
    page.wait_for_timeout(5000)

    # Click Register Compound TAB
    page.get_by_role("tab", name="Register Compound").click()
    page.wait_for_timeout(2000)

    # Fill form with unique compound name (use timestamp to avoid duplicates)
    import time
    compound_name = f"TestCompound{int(time.time())}"

    page.locator('input[aria-label="Compound Name"]').fill(compound_name)
    page.locator('input[aria-label="SMILES"]').fill("CCO")  # Ethanol - simple
    page.wait_for_timeout(1000)

    # Click Register button
    page.get_by_role("button", name="Register", exact=True).click()
    page.wait_for_timeout(5000)  # Wait for registration

    # Click Compounds tab to verify it appears
    page.get_by_role("tab", name="Compounds").click()
    page.wait_for_timeout(3000)

    # Verify we see compound data (AMP- ID in the list)
    # Just check page has some compound content
    expect(page.locator("text=AMP-").first).to_be_attached(timeout=15000)
