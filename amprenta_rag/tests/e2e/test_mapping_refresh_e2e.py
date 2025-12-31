"""End-to-end tests for Mapping Refresh page."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


class TestMappingRefreshE2E:
    """Test Mapping Refresh page functionality."""
    
    @pytest.fixture(autouse=True)
    def setup_page(self, page: Page, streamlit_server: str):
        """Navigate to Mapping Refresh page before each test."""
        page.goto(f"{streamlit_server}/?page=Mapping+Refresh")
        page.wait_for_load_state("domcontentloaded")
        
        # Wait for page to load
        page.wait_for_selector("[data-testid='stMainBlockContainer']", timeout=10000)
    
    def test_page_loads(self, page: Page):
        """Test that page renders with header and 4 tabs."""
        # Check header
        expect(page.get_by_text("🔄 ID Mapping Refresh")).to_be_visible()
        expect(page.get_by_text("Monitor UniProt/KEGG ID mappings")).to_be_visible()
        
        # Check tabs
        expect(page.get_by_role("tab", name="Status")).to_be_visible()
        expect(page.get_by_role("tab", name="Statistics")).to_be_visible()
        expect(page.get_by_role("tab", name="Lookup")).to_be_visible()
        expect(page.get_by_role("tab", name="Jobs")).to_be_visible()
    
    def test_status_tab_shows_metrics(self, page: Page):
        """Test that status tab loads and shows mapping metrics."""
        # Click on Status tab (should be default)
        page.get_by_role("tab", name="Status").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Status tab content
        status_tab = page.get_by_label("Status")
        
        # Check for status subheader
        expect(status_tab.get_by_text("📊 Mapping Status")).to_be_visible()
        
        # Check for metric labels (these should appear even if API is down)
        expect(status_tab.get_by_text("Total Mappings")).to_be_visible()
        expect(status_tab.get_by_text("Expired")).to_be_visible()
        expect(status_tab.get_by_text("Types")).to_be_visible()
        
        # Check for admin section text (since mock user is admin, should see admin controls)
        expect(status_tab.get_by_text("Manual Refresh (Admin Only)")).to_be_visible()
    
    def test_statistics_tab_shows_charts(self, page: Page):
        """Test that statistics tab renders chart sections."""
        # Click on Statistics tab
        page.get_by_role("tab", name="Statistics").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Statistics tab content
        stats_tab = page.get_by_label("Statistics")
        
        # Check for statistics subheader
        expect(stats_tab.get_by_text("📈 Detailed Statistics")).to_be_visible()
        
        # Check for chart section headers
        expect(stats_tab.get_by_text("📊 Coverage by Source Type")).to_be_visible()
        expect(stats_tab.get_by_text("🎯 Coverage by Target Type")).to_be_visible()
        
        # Check for metrics
        expect(stats_tab.get_by_text("Total Mappings")).to_be_visible()
        expect(stats_tab.get_by_text("Permanent Mappings")).to_be_visible()
        expect(stats_tab.get_by_text("TTL Mappings")).to_be_visible()
    
    def test_lookup_single_id(self, page: Page):
        """Test single ID lookup functionality."""
        # Click on Lookup tab
        page.get_by_role("tab", name="Lookup").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Lookup tab content
        lookup_tab = page.get_by_label("Lookup")
        
        # Check for lookup subheader
        expect(lookup_tab.get_by_text("🔍 ID Lookup")).to_be_visible()
        expect(lookup_tab.get_by_text("Single ID Lookup")).to_be_visible()
        
        # Find form elements within the tab
        source_type_select = lookup_tab.locator("[data-testid='stSelectbox']").first
        expect(source_type_select).to_be_visible()
        
        source_id_input = lookup_tab.get_by_placeholder("e.g., TP53, BRCA1")
        expect(source_id_input).to_be_visible()
        
        # Fill in the form
        source_id_input.fill("TP53")
        page.keyboard.press("Tab")
        page.wait_for_timeout(1000)
        
        # Check fallback checkbox
        fallback_checkbox = lookup_tab.get_by_text("Enable API fallback")
        expect(fallback_checkbox).to_be_visible()
        
        # Find and check the lookup button
        lookup_button = lookup_tab.get_by_role("button", name="🔍 Look Up")
        expect(lookup_button).to_be_visible()
        expect(lookup_button).to_be_enabled()
    
    def test_batch_lookup_validation(self, page: Page):
        """Test batch lookup with limit enforcement."""
        # Click on Lookup tab
        page.get_by_role("tab", name="Lookup").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Lookup tab content
        lookup_tab = page.get_by_label("Lookup")
        
        # Check for batch lookup section
        expect(lookup_tab.get_by_text("Batch Lookup")).to_be_visible()
        expect(lookup_tab.get_by_text("maximum 1000 IDs per request")).to_be_visible()
        
        # Find batch form elements
        batch_source_select = lookup_tab.locator("[data-testid='stSelectbox']").nth(1)
        expect(batch_source_select).to_be_visible()
        
        batch_target_select = lookup_tab.locator("[data-testid='stSelectbox']").nth(2)
        expect(batch_target_select).to_be_visible()
        
        ids_textarea = lookup_tab.get_by_placeholder("TP53\nBRCA1\nEGFR")
        expect(ids_textarea).to_be_visible()
        
        # Fill in some test IDs
        ids_textarea.fill("TP53\nBRCA1\nEGFR")
        page.keyboard.press("Tab")
        page.wait_for_timeout(1000)
        
        # Find batch lookup button
        batch_button = lookup_tab.get_by_role("button", name="🔍 Batch Look Up")
        expect(batch_button).to_be_visible()
        expect(batch_button).to_be_enabled()
    
    def test_jobs_tab_shows_content(self, page: Page):
        """Test jobs tab displays job information."""
        # Click on Jobs tab
        page.get_by_role("tab", name="Jobs").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Jobs tab content
        jobs_tab = page.get_by_label("Jobs")
        
        # Check for jobs subheader
        expect(jobs_tab.get_by_text("🔄 Refresh Jobs")).to_be_visible()
        
        # Check for job history info
        expect(jobs_tab.get_by_text("Job history tracking coming soon")).to_be_visible()
        
        # Check for current status section
        expect(jobs_tab.get_by_text("Current Status")).to_be_visible()
        expect(jobs_tab.get_by_text("Scheduled Jobs")).to_be_visible()
        expect(jobs_tab.get_by_text("UniProt mappings refresh weekly")).to_be_visible()
        
        # Check for admin quick actions (since mock user is admin)
        expect(jobs_tab.get_by_text("Quick Actions (Admin)")).to_be_visible()
    
    def test_refresh_requires_admin(self, page: Page):
        """Test that refresh functionality shows admin controls for admin user."""
        # Click on Status tab
        page.get_by_role("tab", name="Status").click()
        page.wait_for_timeout(1000)
        
        # Scope to the Status tab content
        status_tab = page.get_by_label("Status")
        
        # Check that admin controls are shown (since mock user is admin)
        expect(status_tab.get_by_text("Manual Refresh (Admin Only)")).to_be_visible()
        expect(status_tab.get_by_role("button", name="🔄 Refresh UniProt Mappings")).to_be_visible()
        
        # Switch to Jobs tab and check admin controls there too
        page.get_by_role("tab", name="Jobs").click()
        page.wait_for_timeout(1000)
        
        jobs_tab = page.get_by_label("Jobs")
        expect(jobs_tab.get_by_text("Quick Actions (Admin)")).to_be_visible()
        expect(jobs_tab.get_by_role("button", name="🔄 Trigger UniProt Refresh")).to_be_visible()
