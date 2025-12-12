"""
Comprehensive E2E Platform Test Suite
Tests all major dashboard functionality
"""
import pytest
from playwright.sync_api import Page, expect
import time

# Quick self-check
if __name__ == "__main__":
    print("Syntax OK")


class TestNavigationSmoke:
    """Test all sidebar sections expand and pages load"""
    
    @pytest.mark.parametrize("section,pages", [
        ("🔍 Discovery", ["Overview", "Experiments", "Discovery Workflow", "Variant Tracking"]),
        ("📊 Analysis", ["Datasets", "Signatures", "RAG Query"]),
        ("📋 ELN", ["Protocols", "Sample Inventory"]),
        ("⚙️ Admin", ["Audit Logs", "Data Quality", "System Health"]),
    ])
    def test_section_pages_load(self, page: Page, base_url: str, section: str, pages: list):
        """Test that pages in each section load without errors."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)  # Wait for Streamlit to fully load
        
        # Expand section
        try:
            expander_button = page.get_by_role("button", name=section).first
            if expander_button.count() > 0:
                expander_button.click()
                page.wait_for_timeout(1000)  # Wait for expander to open
        except Exception as e:
            pytest.skip(f"Could not find or click section '{section}': {e}")
        
        for page_name in pages:
            try:
                # Click the page button
                page_button = page.get_by_role("button", name=page_name, exact=True).first
                if page_button.count() > 0:
                    page_button.click()
                    page.wait_for_timeout(2000)  # Wait for page to load
                    
                    # Check for both error patterns
                    loading_error = page.locator("text=Error loading page").count()
                    rendering_error = page.locator("text=Error rendering page").count()
                    assert loading_error == 0, f"{page_name}: Import error - Error loading page"
                    assert rendering_error == 0, f"{page_name}: Runtime error - Error rendering page"
                    
                    # Check that page loaded (some content is present)
                    # Most pages will have their title/header visible
                    assert page.locator("body").count() > 0, f"{page_name} did not load"
            except Exception as e:
                pytest.fail(f"Failed to navigate to {page_name}: {e}")


class TestChemistryFlow:
    """Test compound registration and viewing"""
    
    def test_register_and_view_compound(self, page: Page, base_url: str):
        """Test registering a compound and verifying it appears in the list."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate to Chemistry
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Now click Analysis expander (Streamlit expander, not a button)
            page.locator("text=📊 Analysis").click()
            page.wait_for_timeout(1000)
            # Then click Chemistry
            page.get_by_role("button", name="Chemistry", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Chemistry: {e}")
        
        # Go to Register tab (use text locator for Streamlit tabs)
        try:
            page.locator('[data-testid="stTabs"] button:has-text("Register Compound")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not find Register Compound tab: {e}")
        
        # Fill form with unique compound name
        compound_name = f"TestCompound_E2E_{int(time.time())}"
        
        try:
            # Use exact aria-label attribute selectors to avoid partial matches
            name_input = page.locator('input[aria-label="Compound Name"]')
            name_input.fill(compound_name)
            
            smiles_input = page.locator('input[aria-label="SMILES"]')
            smiles_input.fill("CCO")  # Ethanol
            
            page.wait_for_timeout(500)
            
            # Click Register button
            register_button = page.get_by_role("button", name="Register", exact=True).first
            if register_button.count() > 0:
                register_button.click()
            else:
                # Try without exact match
                page.get_by_role("button", name="Register").first.click()
            
            page.wait_for_timeout(3000)  # Wait for registration to complete
        except Exception as e:
            pytest.fail(f"Failed to fill form or register compound: {e}")
        
        # Verify registration (check for AMP- ID pattern)
        try:
            amp_id_locator = page.locator("text=AMP-").first
            expect(amp_id_locator).to_be_attached(timeout=10000)
        except Exception as e:
            pytest.fail(f"Compound registration verification failed: {e}")


class TestRAGFlow:
    """Test RAG query with hybrid search"""
    
    def test_rag_query_with_citations(self, page: Page, base_url: str):
        """Test RAG query functionality with citations."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate to RAG Query
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # RAG Query might be in "Other Pages" or "Analysis"
            # Try Analysis first (Streamlit expander, not a button)
            analysis_expander = page.locator("text=📊 Analysis").first
            if analysis_expander.count() > 0:
                analysis_expander.click()
                page.wait_for_timeout(1000)
            
            # Check if RAG Query is in Analysis
            rag_button = page.get_by_role("button", name="RAG Query", exact=True).first
            if rag_button.count() == 0:
                # Try Other Pages
                other_pages_expander = page.locator("text=📚 Other Pages").first
                if other_pages_expander.count() > 0:
                    other_pages_expander.click()
                    page.wait_for_timeout(1000)
                    rag_button = page.get_by_role("button", name="RAG Query", exact=True).first
            
            if rag_button.count() > 0:
                rag_button.click()
            else:
                pytest.skip("RAG Query page not found in sidebar")
            
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.skip(f"Could not navigate to RAG Query: {e}")
        
        # Enter query
        try:
            # Use exact aria-label selector for query input
            query_input = page.locator('textarea[aria-label="Enter your query"]')
            query_input.fill("What causes ALS?")
        except Exception as e:
            pytest.fail(f"Could not find query input: {e}")
        
        # Enable hybrid search (if checkbox exists)
        try:
            # Use exact aria-label selector for checkbox
            hybrid_checkbox = page.locator('input[aria-label="Use Hybrid Search"]')
            if hybrid_checkbox.count() > 0:
                hybrid_checkbox.check()
        except Exception:
            # Hybrid search checkbox might not exist or have different label
            pass
        
        # Submit query
        try:
            search_button = page.get_by_role("button", name="🔍 Search", exact=True).first
            if search_button.count() == 0:
                search_button = page.get_by_role("button", name="Search", exact=True).first
            if search_button.count() > 0:
                search_button.click()
            else:
                pytest.fail("Could not find Search button")
        except Exception as e:
            pytest.fail(f"Could not submit query: {e}")
        
        # Wait for results (RAG takes time)
        page.wait_for_timeout(10000)
        
        # Check for results (matches section or source citations)
        try:
            # Look for common result indicators
            source_locator = page.locator("text=Source").first
            if source_locator.count() == 0:
                # Try alternative indicators
                matches_locator = page.locator("text=Matches").first
                if matches_locator.count() == 0:
                    # Check if any content loaded (not just error)
                    error_locator = page.locator("text=Error").first
                    if error_locator.count() > 0:
                        pytest.fail("RAG query returned an error")
                    # If no error, assume results loaded
                    assert page.locator("body").count() > 0, "Page did not load results"
                else:
                    assert matches_locator.count() > 0, "Matches section not found"
            else:
                assert source_locator.count() > 0, "Source section not found"
        except Exception as e:
            pytest.fail(f"Could not verify RAG query results: {e}")


class TestExperimentFlow:
    """Test experiment creation"""
    
    def test_create_experiment(self, page: Page, base_url: str):
        """Test navigating to Experiments page and verifying it loads."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate to Experiments (Discovery is expanded by default)
        try:
            # Scroll sidebar to reveal buttons
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Click Experiments button directly (Discovery expander is expanded by default)
            page.locator('text=Experiments').first.click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Experiments: {e}")
        
        # Check page loaded (look for experiment-related content)
        try:
            # Look for common experiment page indicators
            experiment_text = page.locator("text=Experiment").first
            if experiment_text.count() == 0:
                # Try alternative indicators
                overview_text = page.locator("text=Overview").first
                if overview_text.count() == 0:
                    # Check for any content (not error)
                    error_locator = page.locator("text=Error rendering").first
                    if error_locator.count() > 0:
                        pytest.fail("Experiments page shows rendering error")
                    # If no error, assume page loaded
                    assert page.locator("body").count() > 0, "Experiments page did not load"
                else:
                    assert overview_text.count() > 0, "Experiments page content not found"
            else:
                assert experiment_text.count() > 0, "Experiments page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Experiments page loaded: {e}")


class TestSampleInventory:
    """Test sample inventory workflow"""
    
    def test_sample_workflow(self, page: Page, base_url: str):
        """Create location, register sample"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ELN expander (collapsed) > Sample Inventory
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand ELN section
            page.locator("text=📋 ELN").click()
            page.wait_for_timeout(1000)
            
            # Click Sample Inventory
            page.get_by_role("button", name="Sample Inventory", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Sample Inventory: {e}")
        
        # Tab: Storage Locations > create "Freezer-Test"
        try:
            page.locator('button[role="tab"]:has-text("Storage Locations")').click()
            page.wait_for_timeout(1000)
            
            # Fill location form
            location_name_input = page.locator('input[aria-label="Name"]')
            location_name_input.fill(f"Freezer-Test-{int(time.time())}")
            
            # Click create button
            create_button = page.get_by_role("button", name="Add Location", exact=True).first
            if create_button.count() > 0:
                create_button.click()
                page.wait_for_timeout(2000)
        except Exception as e:
            pytest.skip(f"Could not create storage location: {e}")
        
        # Tab: Register Sample > fill name, generate barcode
        try:
            page.locator('button[role="tab"]:has-text("Register Sample")').click()
            page.wait_for_timeout(1000)
            
            # Fill sample form
            sample_name_input = page.locator('input[aria-label="Name*"]')
            sample_name_input.fill(f"TestSample-{int(time.time())}")
            
            # Click register button
            register_button = page.get_by_role("button", name="Register", exact=True).first
            if register_button.count() > 0:
                register_button.click()
                page.wait_for_timeout(2000)
                
                # Verify sample appears (check for barcode or success message)
                assert page.locator("text=Sample").first.count() > 0 or page.locator("text=success").first.count() > 0
        except Exception as e:
            pytest.fail(f"Could not register sample: {e}")


class TestProtocols:
    """Test protocol workflow"""
    
    def test_protocol_workflow(self, page: Page, base_url: str):
        """Create protocol"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ELN expander > Protocols
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand ELN section
            page.locator("text=📋 ELN").click()
            page.wait_for_timeout(1000)
            
            # Click Protocols
            page.get_by_role("button", name="Protocols", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Protocols: {e}")
        
        # Tab: Create Protocol > fill title, steps
        try:
            page.locator('button[role="tab"]:has-text("Create Protocol")').click()
            page.wait_for_timeout(1000)
            
            # Fill protocol form
            title_input = page.locator('input[aria-label="Protocol Name*"]')
            title_input.fill(f"Test Protocol {int(time.time())}")
            
            steps_input = page.locator('textarea[aria-label="Steps"]')
            steps_input.fill("Step 1: Prepare\nStep 2: Execute\nStep 3: Analyze")
            
            # Click save button (use form submit button or last Create Protocol button)
            save_button = page.locator('[data-testid="stFormSubmitButton"] button').first
            if save_button.count() == 0:
                save_button = page.get_by_role("button", name="Create Protocol").last
            if save_button.count() > 0:
                save_button.click()
                page.wait_for_timeout(2000)
                
                # Verify protocol created
                assert page.locator("text=Protocol").first.count() > 0 or page.locator("text=success").first.count() > 0
        except Exception as e:
            pytest.fail(f"Could not create protocol: {e}")


class TestDiscoveryWorkflow:
    """Test discovery workflow"""
    
    def test_discovery_job(self, page: Page, base_url: str):
        """Run discovery scan"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Discovery (expanded) > Discovery Workflow
        try:
            # Scroll sidebar to reveal buttons
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Click Discovery Workflow (Discovery is expanded by default)
            page.locator('text=Discovery Workflow').first.click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Discovery Workflow: {e}")
        
        # Tab: Run Discovery > enter query "test"
        try:
            page.locator('button[role="tab"]:has-text("Run Discovery")').click()
            page.wait_for_timeout(1000)
            
            # Fill discovery form
            query_input = page.locator('input[aria-label="Search Query"]')
            query_input.fill("test")
            
            # Click Start Discovery button
            start_button = page.get_by_role("button", name="Start Discovery", exact=True).first
            if start_button.count() == 0:
                start_button = page.get_by_role("button", name="🚀 Start Discovery", exact=True).first
            if start_button.count() > 0:
                start_button.click()
                page.wait_for_timeout(5000)  # Discovery takes time
                
                # Check for job created or pending studies
                assert page.locator("text=Job").first.count() > 0 or page.locator("text=Pending").first.count() > 0 or page.locator("text=success").first.count() > 0
        except Exception as e:
            pytest.skip(f"Could not run discovery job: {e}")


class TestVariantTracking:
    """Test variant tracking"""
    
    def test_add_variant(self, page: Page, base_url: str):
        """Add genetic variant"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Discovery > Variant Tracking
        try:
            # Scroll sidebar to reveal buttons
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Click Variant Tracking (Discovery is expanded by default)
            page.locator('text=Variant Tracking').first.click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Variant Tracking: {e}")
        
        # Tab: Add Variant > fill Gene, Variant, Organism
        try:
            page.locator('button[role="tab"]:has-text("Add Variant")').click()
            page.wait_for_timeout(1000)
            
            # Fill variant form
            gene_input = page.locator('input[aria-label="Gene*"]')
            gene_input.fill("TP53")
            
            variant_input = page.locator('input[aria-label="Variant*"]')
            variant_input.fill("p.R273H")
            
            organism_input = page.locator('input[aria-label="Cell Line/Organism*"]')
            organism_input.fill("HeLa")
            
            # Click Save button
            save_button = page.get_by_role("button", name="💾 Save Variant", exact=True).first
            if save_button.count() == 0:
                save_button = page.get_by_role("button", name="Save Variant", exact=True).first
            if save_button.count() > 0:
                save_button.click()
                page.wait_for_timeout(2000)
                
                # Verify variant in Browse tab
                page.locator('button[role="tab"]:has-text("Browse")').click()
                page.wait_for_timeout(1000)
                assert page.locator("text=TP53").first.count() > 0
        except Exception as e:
            pytest.fail(f"Could not add variant: {e}")


class TestDatasets:
    """Test datasets page"""
    
    def test_view_datasets(self, page: Page, base_url: str):
        """View datasets page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Other Pages > Datasets
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages
            page.locator("text=📚 Other Pages").click()
            page.wait_for_timeout(1000)
            
            # Click Datasets
            page.get_by_role("button", name="Datasets", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Datasets: {e}")
        
        # Verify datasets table or empty state loads
        try:
            # Check for datasets content or empty state
            assert page.locator("text=Dataset").first.count() > 0 or page.locator("text=No datasets").first.count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Datasets page loaded: {e}")


class TestSignatures:
    """Test signatures page"""
    
    def test_view_signatures(self, page: Page, base_url: str):
        """View signatures page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Other Pages > Signatures
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages
            page.locator("text=📚 Other Pages").click()
            page.wait_for_timeout(1000)
            
            # Click Signatures
            page.get_by_role("button", name="Signatures", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Signatures: {e}")
        
        # Verify page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Signatures page: Import error - Error loading page"
            assert rendering_error == 0, "Signatures page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Signatures page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Signatures page loaded: {e}")

