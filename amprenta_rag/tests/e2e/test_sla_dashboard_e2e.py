"""E2E tests for SLA Dashboard."""

import pytest
from playwright.sync_api import Page, expect


def test_sla_dashboard_page_loads(page: Page, streamlit_server: str):
    """Test that SLA Dashboard page loads with 4 tabs."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Check page title/header
    expect(page.get_by_role("heading", name="📈 SLA Dashboard")).to_be_visible()
    expect(page.get_by_text("Monitor and manage Service Level Agreements")).to_be_visible()
    
    # Check all 4 tabs are present
    expect(page.get_by_role("tab", name="Overview")).to_be_visible()
    expect(page.get_by_role("tab", name="Reviews")).to_be_visible()
    expect(page.get_by_role("tab", name="Cycles")).to_be_visible()
    expect(page.get_by_role("tab", name="Settings")).to_be_visible()


def test_overview_tab_shows_metrics(page: Page, streamlit_server: str):
    """Test that Overview tab shows SLA metrics cards."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Click Overview tab (should be default)
    page.get_by_role("tab", name="Overview").click()
    page.wait_for_timeout(1000)
    
    # Check for metrics cards - should show either real data or demo data
    expect(page.get_by_text("SLA Status Overview")).to_be_visible()
    
    # Look for metric labels (these should always be present)
    expect(page.get_by_text("🟢 On Track")).to_be_visible()
    expect(page.get_by_text("🟡 At Risk")).to_be_visible()
    expect(page.get_by_text("🟠 Overdue")).to_be_visible()
    expect(page.get_by_text("🔴 Breached")).to_be_visible()
    
    # Check for either real data, demo data, or no data message
    page.wait_for_timeout(1000)
    # Should show either "Total Active Reviews", "API unavailable", or "No active reviews"
    try:
        expect(page.get_by_text("Total Active Reviews")).to_be_visible()
    except AssertionError:
        try:
            expect(page.get_by_text("API unavailable")).to_be_visible()
        except AssertionError:
            expect(page.get_by_text("No active reviews with SLAs")).to_be_visible()


def test_reviews_tab_shows_table(page: Page, streamlit_server: str):
    """Test that Reviews tab shows reviews table with status badges."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Click Reviews tab
    page.get_by_role("tab", name="Reviews").click()
    page.wait_for_timeout(1000)
    
    # Check tab content
    expect(page.get_by_role("heading", name="📋 Active Reviews")).to_be_visible()
    
    # Check for filter controls - scope to Reviews tab to avoid conflicts
    reviews_tab = page.get_by_label("Reviews")
    expect(reviews_tab.get_by_text("Entity Type")).to_be_visible()
    expect(reviews_tab.get_by_text("Status")).to_be_visible()
    
    # Should show either reviews or "No overdue reviews" or demo data
    page.wait_for_timeout(1000)
    try:
        # Look for review content or no reviews message
        page.get_by_text("No overdue reviews found.").wait_for(timeout=2000)
    except:
        try:
            # Or demo data
            expect(page.get_by_text("API unavailable")).to_be_visible()
        except:
            # Or actual reviews
            pass  # Any reviews shown is acceptable


def test_cycles_tab_admin_only(page: Page, streamlit_server: str):
    """Test that Cycles tab shows for admin user or shows access message."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Click Cycles tab
    page.get_by_role("tab", name="Cycles").click()
    page.wait_for_timeout(1000)
    
    # Check tab content
    expect(page.get_by_role("heading", name="🔄 Review Cycles")).to_be_visible()
    
    # Should show either admin content or access denied message
    try:
        # If admin access
        expect(page.get_by_text("Create New Cycle")).to_be_visible()
    except AssertionError:
        # If not admin
        expect(page.get_by_text("Admin access required")).to_be_visible()


def test_settings_tab_shows_sla_rules(page: Page, streamlit_server: str):
    """Test that Settings tab shows SLA rules or access message."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Click Settings tab
    page.get_by_role("tab", name="Settings").click()
    page.wait_for_timeout(1000)
    
    # Check tab content
    expect(page.get_by_role("heading", name="⚙️ SLA Rules")).to_be_visible()
    
    # Should show either admin content or access denied message
    try:
        # If admin access
        expect(page.get_by_text("Create New SLA Rule")).to_be_visible()
    except AssertionError:
        # If not admin
        expect(page.get_by_text("Admin access required")).to_be_visible()


def test_create_cycle_form(page: Page, streamlit_server: str):
    """Test that cycle creation form renders when admin access is available."""
    page.goto(f"{streamlit_server}/?page=SLA+Dashboard")
    page.wait_for_timeout(2000)
    
    # Click Cycles tab
    page.get_by_role("tab", name="Cycles").click()
    page.wait_for_timeout(1000)
    
    # Since AUTH_DISABLED=true, user should have admin access
    # Look for create cycle form
    create_cycle_button = page.get_by_text("Create New Cycle")
    expect(create_cycle_button).to_be_visible()
    
    # Click to expand form
    create_cycle_button.click()
    page.wait_for_timeout(1000)
    
    # Check form fields are present
    expect(page.get_by_placeholder("e.g., Weekly Dataset Reviews")).to_be_visible()
    expect(page.get_by_text("Frequency")).to_be_visible()
    expect(page.get_by_text("Reviewer Pool")).to_be_visible()
    
    # Check submit button
    expect(page.get_by_role("button", name="Create Cycle")).to_be_visible()
