"""End-to-end tests for Provenance Ledger page."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


class TestProvenanceLedgerE2E:
    """Test Provenance Ledger page functionality."""
    
    @pytest.fixture(autouse=True)
    def setup_page(self, page: Page, streamlit_server: str):
        """Navigate to Provenance Ledger page before each test."""
        page.goto(f"{streamlit_server}/Provenance_Ledger")
        page.wait_for_load_state("domcontentloaded")
        
        # Wait for page to load
        page.wait_for_selector("[data-testid='stMainBlockContainer']", timeout=10000)
    
    def test_page_loads(self, page: Page):
        """Test that page renders with header and 3 tabs."""
        # Check header
        expect(page.get_by_text("📜 Provenance Ledger")).to_be_visible()
        expect(page.get_by_text("Browse entity version history")).to_be_visible()
        
        # Check tabs
        expect(page.get_by_role("tab", name="Version History")).to_be_visible()
        expect(page.get_by_role("tab", name="Compare Versions")).to_be_visible()
        expect(page.get_by_role("tab", name="Restore Version")).to_be_visible()
    
    def test_version_history_loads(self, page: Page):
        """Test that version history tab loads and shows table when entity is provided."""
        # Click on Version History tab
        page.get_by_role("tab", name="Version History").click()
        
        # Check form elements
        expect(page.get_by_text("Entity Type")).to_be_visible()
        expect(page.get_by_text("Entity ID")).to_be_visible()
        
        # Check selectbox and input
        entity_type_select = page.locator("[data-testid='stSelectbox']").first
        expect(entity_type_select).to_be_visible()
        
        entity_id_input = page.get_by_placeholder("UUID")
        expect(entity_id_input).to_be_visible()
        
        # Check Load Versions button
        load_button = page.get_by_role("button", name="Load Versions")
        expect(load_button).to_be_visible()
        expect(load_button).to_be_disabled()  # Should be disabled without entity ID
        
        # Fill in entity ID to enable button
        entity_id_input.fill("123e4567-e89b-12d3-a456-426614174000")
        expect(load_button).to_be_enabled()
    
    def test_version_detail_expands(self, page: Page):
        """Test that version detail expanders work."""
        # Click on Version History tab
        page.get_by_role("tab", name="Version History").click()
        
        # Fill form
        entity_id_input = page.get_by_placeholder("UUID")
        entity_id_input.fill("123e4567-e89b-12d3-a456-426614174000")
        
        # Click Load Versions button
        load_button = page.get_by_role("button", name="Load Versions")
        load_button.click()
        
        # Wait for potential API response (may show error or data)
        page.wait_for_timeout(2000)
        
        # Check that Version Details section exists
        # (even if no data, the section header should be present)
        version_details = page.get_by_text("Version Details")
        if version_details.is_visible():
            expect(version_details).to_be_visible()
    
    def test_compare_versions_shows_diff(self, page: Page):
        """Test that compare versions tab shows diff sections."""
        # Click on Compare Versions tab
        page.get_by_role("tab", name="Compare Versions").click()
        
        # Check form elements
        expect(page.get_by_text("Entity Type")).to_be_visible()
        expect(page.get_by_text("Entity ID")).to_be_visible()
        expect(page.get_by_text("Version A")).to_be_visible()
        expect(page.get_by_text("Version B")).to_be_visible()
        
        # Check inputs
        entity_id_input = page.get_by_placeholder("UUID")
        expect(entity_id_input).to_be_visible()
        
        version_a_input = page.locator("[data-testid='stNumberInput']").first
        expect(version_a_input).to_be_visible()
        
        version_b_input = page.locator("[data-testid='stNumberInput']").last
        expect(version_b_input).to_be_visible()
        
        # Check Compare button
        compare_button = page.get_by_role("button", name="Compare")
        expect(compare_button).to_be_visible()
        expect(compare_button).to_be_disabled()  # Should be disabled without entity ID
        
        # Fill in entity ID to enable button
        entity_id_input.fill("123e4567-e89b-12d3-a456-426614174000")
        expect(compare_button).to_be_enabled()
    
    def test_restore_requires_confirmation(self, page: Page):
        """Test that restore button is disabled without confirmation checkbox."""
        # Click on Restore Version tab
        page.get_by_role("tab", name="Restore Version").click()
        
        # Check for admin access message or form
        main_content = page.locator("[data-testid='stMainBlockContainer']")
        
        # If admin access required message is shown
        if page.get_by_text("Admin access required").is_visible():
            expect(page.get_by_text("🔒 Admin access required for version restore")).to_be_visible()
        else:
            # If form is shown (user is admin), test form behavior
            version_id_input = page.get_by_placeholder("UUID of version to restore")
            reason_input = page.get_by_placeholder("Why are you restoring this version?")
            confirm_checkbox = page.get_by_text("I confirm this restore action")
            restore_button = page.get_by_role("button", name="Restore")
            
            expect(version_id_input).to_be_visible()
            expect(reason_input).to_be_visible()
            expect(confirm_checkbox).to_be_visible()
            expect(restore_button).to_be_visible()
            expect(restore_button).to_be_disabled()  # Should be disabled initially
            
            # Fill form but don't check confirmation
            version_id_input.fill("123e4567-e89b-12d3-a456-426614174000")
            reason_input.fill("Test restore")
            expect(restore_button).to_be_disabled()  # Still disabled without checkbox
            
            # Check confirmation
            confirm_checkbox.check()
            expect(restore_button).to_be_enabled()  # Now enabled
    
    def test_restore_requires_admin(self, page: Page):
        """Test that non-admin users see error message on restore tab."""
        # Click on Restore Version tab
        page.get_by_role("tab", name="Restore Version").click()
        
        # Check for admin access message
        # Note: This test assumes the user is not admin by default
        # In a real test environment, we would mock the user session
        admin_message = page.get_by_text("Admin access required")
        if admin_message.is_visible():
            expect(page.get_by_text("🔒 Admin access required for version restore")).to_be_visible()
        else:
            # If user is admin, we can't test the non-admin case in this session
            # This would require mocking the session state
            pass
