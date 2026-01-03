"""E2E tests for Report Builder dashboard."""

import pytest
from playwright.sync_api import Page, expect


@pytest.fixture
def report_builder_page(page: Page, streamlit_server: str):
    """Navigate to Report Builder page."""
    # Navigate directly to Report Builder page via URL parameter
    page.goto(f"{streamlit_server}/?page=Report+Builder")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)  # Allow page to fully render
    return page


def test_page_loads(report_builder_page: Page):
    """Report Builder page loads with main title."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_role("heading", name="📄 Custom Report Builder")).to_be_visible()


def test_build_report_tab_loads(report_builder_page: Page):
    """Build Report tab loads with section library."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Should be on Build Report tab by default
    expect(main.get_by_text("➕ Add Sections")).to_be_visible()
    expect(main.get_by_text("📋 Report Sections")).to_be_visible()
    expect(main.get_by_text("Add sections from the left panel to build your report.")).to_be_visible()


def test_templates_tab_loads(report_builder_page: Page):
    """Templates tab loads with save button."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click My Templates tab
    main.get_by_role("tab", name="📁 My Templates").click()
    report_builder_page.wait_for_timeout(2000)
    
    expect(main.get_by_text("💾 Save Current as Template")).to_be_visible()
    expect(main.get_by_text("📁 Saved Templates")).to_be_visible()


def test_preview_tab_loads(report_builder_page: Page):
    """Preview tab loads with generate button."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Preview tab
    main.get_by_role("tab", name="👁️ Preview").click()
    report_builder_page.wait_for_timeout(2000)
    
    expect(main.get_by_text("🔄 Generate Preview")).to_be_visible()
    expect(main.get_by_text("👁️ Report Preview")).to_be_visible()


def test_export_tab_loads(report_builder_page: Page):
    """Export tab loads with controls."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Export tab
    main.get_by_role("tab", name="📥 Export").click()
    report_builder_page.wait_for_timeout(2000)
    
    expect(main.get_by_text("📥 Export Report")).to_be_visible()
    expect(main.get_by_text("Report Title")).to_be_visible()
    expect(main.get_by_text("Format")).to_be_visible()
    expect(main.get_by_text("📥 Generate & Download")).to_be_visible()


def test_section_library_tab_loads(report_builder_page: Page):
    """Section Library tab shows all section types."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Click Section Library tab
    main.get_by_role("tab", name="📚 Section Library").click()
    report_builder_page.wait_for_timeout(2000)
    
    expect(main.get_by_text("📚 Available Section Types")).to_be_visible()
    expect(main.get_by_text("Browse all available section types")).to_be_visible()


def test_navigation_between_tabs(report_builder_page: Page):
    """Test navigation between all tabs works correctly."""
    main = report_builder_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Test all tab clicks
    tabs = ["📁 My Templates", "👁️ Preview", "📥 Export", "📚 Section Library", "🔨 Build Report"]
    
    for tab_name in tabs:
        main.get_by_role("tab", name=tab_name).click()
        report_builder_page.wait_for_timeout(1000)
        # Just verify we can click each tab without errors
        expect(main).to_be_visible()
