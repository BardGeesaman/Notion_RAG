"""E2E tests for Compound Inventory dashboard."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def inventory_page(page: Page, streamlit_server: str):
    """Navigate to Compound Inventory page."""
    # Navigate directly to Compound Inventory page via URL parameter
    page.goto(f"{streamlit_server}/?page=Compound+Inventory")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)  # Allow page to fully render
    return page


def test_compound_stocks_tab_loads(inventory_page: Page):
    """Compound Stocks tab loads with table and filters."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Check main title
    expect(main.get_by_role("heading", name="🧪 Compound Inventory")).to_be_visible()
    
    # Check tab content (Compound Stocks tab should be active by default)
    expect(main.get_by_role("heading", name="Compound Sample Inventory")).to_be_visible()
    # Use more specific selector to avoid strict mode violation
    expect(main.get_by_text("Low Stock Only")).to_be_visible()


def test_plates_tab_loads(inventory_page: Page):
    """Plates tab loads with table and New Plate button."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Plates tab
    main.get_by_role("tab", name="🔬 Plates").click()
    inventory_page.wait_for_timeout(2000)
    
    # Check plates content
    expect(main.get_by_role("heading", name="Compound Plates")).to_be_visible()
    expect(main.get_by_text("New Plate")).to_be_visible()


def test_requests_tab_loads(inventory_page: Page):
    """Requests tab loads with view modes and New Request button."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Requests tab
    main.get_by_role("tab", name="📋 Requests").click()
    inventory_page.wait_for_timeout(2000)
    
    # Check requests content
    expect(main.get_by_role("heading", name="Compound Requests")).to_be_visible()
    expect(main.get_by_text("New Request")).to_be_visible()


def test_register_sample_form_visible(inventory_page: Page):
    """Register Sample tab has form with required fields."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Register Sample tab
    main.get_by_role("tab", name="➕ Register Sample").click()
    inventory_page.wait_for_timeout(2000)
    
    # Check form elements - look for subheaders that should be visible
    expect(main.get_by_text("Compound & Quantity")).to_be_visible()
    expect(main.get_by_text("Storage & Tracking")).to_be_visible()


def test_barcode_lookup_input(inventory_page: Page):
    """Barcode Lookup has input field and recent samples."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Barcode Lookup tab
    main.get_by_role("tab", name="🔍 Barcode Lookup").click()
    inventory_page.wait_for_timeout(2000)
    
    # Check barcode lookup elements
    expect(main.get_by_role("heading", name="🔍 Barcode Lookup")).to_be_visible()
    expect(main.get_by_text("Recent Samples")).to_be_visible()
    expect(main.get_by_placeholder("COMP-20250103-ABC123")).to_be_visible()


def test_navigation_between_tabs(inventory_page: Page):
    """Test navigation between all tabs works correctly."""
    main = inventory_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Test all tab clicks
    tabs = ["🔬 Plates", "📋 Requests", "➕ Register Sample", "🔍 Barcode Lookup", "📦 Compound Stocks"]
    
    for tab_name in tabs:
        main.get_by_role("tab", name=tab_name).click()
        inventory_page.wait_for_timeout(1000)
        # Just verify we can click each tab without errors
        expect(main).to_be_visible()
