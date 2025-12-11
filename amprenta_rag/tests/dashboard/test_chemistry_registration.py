"""Chemistry registration E2E test."""
from playwright.sync_api import Page, expect


def test_register_compound_and_verify(page: Page, base_url: str):
    page.goto(base_url)
    page.wait_for_timeout(5000)
    page.screenshot(path="/tmp/step1_initial.png")
    
    # Click Chemistry in sidebar using radio button
    page.locator("text=Chemistry").click()
    page.wait_for_timeout(5000)
    page.screenshot(path="/tmp/step2_chemistry.png")
    
    # Debug: print page content
    print("Page title:", page.title())
    print("URL:", page.url)
    
    # Try clicking Register Compound tab
    try:
        page.locator("text=Register Compound").click(timeout=10000)
        page.screenshot(path="/tmp/step3_register_tab.png")
    except Exception as e:
        page.screenshot(path="/tmp/step3_error.png")
        raise e
    
    page.wait_for_timeout(2000)
    
    # Fill form
    page.locator('input[aria-label="Compound Name"]').fill("Caffeine")
    page.locator('input[aria-label="SMILES"]').fill("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    page.screenshot(path="/tmp/step4_filled.png")
    
    # Click Register
    page.locator("button").filter(has_text="Register").click()
    page.wait_for_timeout(3000)
    page.screenshot(path="/tmp/step5_registered.png")
    
    # Verify AMP- ID
    expect(page.locator("text=AMP-")).to_be_visible(timeout=10000)
