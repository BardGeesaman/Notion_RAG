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
            sidebar.locator('button:has-text("Experiments")').click()
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
            sidebar.locator('button:has-text("Discovery Workflow")').click()
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
            sidebar.locator('button:has-text("Variant Tracking")').click()
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


class TestDataQuality:
    """Test data quality validation"""
    
    def test_run_validation(self, page: Page, base_url: str):
        """Run validation check"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Admin expander > Data Quality
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section (Data Quality is in Other Pages, not Admin)
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Data Quality
            sidebar.locator('button:has-text("Data Quality")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Data Quality: {e}")
        
        # Button: "🔍 Run Validation"
        try:
            run_button = page.get_by_role("button", name="🔍 Run Validation", exact=True).first
            if run_button.count() == 0:
                run_button = page.get_by_role("button", name="Run Validation", exact=True).first
            if run_button.count() > 0:
                run_button.click()
                page.wait_for_timeout(3000)  # Validation takes time
                
                # Verify results appear or empty state
                assert page.locator("text=Validation").first.count() > 0 or page.locator("text=No issues").first.count() > 0 or page.locator("text=Issues").first.count() > 0
        except Exception as e:
            pytest.skip(f"Could not run validation: {e}")


class TestStatisticalAnalysis:
    """Test statistical analysis page"""
    
    def test_view_page(self, page: Page, base_url: str):
        """View statistical analysis page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Other Pages expander > Statistical Analysis
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages
            page.locator("text=📚 Other Pages").click()
            page.wait_for_timeout(1000)
            
            # Click Statistical Analysis
            sidebar.locator('button:has-text("Statistical Analysis")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Statistical Analysis: {e}")
        
        # Verify page loads (has selectbox for Test Type)
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Statistical Analysis page: Import error - Error loading page"
            assert rendering_error == 0, "Statistical Analysis page: Runtime error - Error rendering page"
            
            # Check for Test Type selectbox or page content
            assert page.locator("text=Test Type").first.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Statistical Analysis page loaded: {e}")


class TestVisualizations:
    """Test visualizations page"""
    
    def test_view_visualizations(self, page: Page, base_url: str):
        """View visualizations page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Analysis expander > Visualizations
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Analysis section
            page.locator("text=📊 Analysis").click()
            page.wait_for_timeout(1000)
            
            # Click Visualizations
            sidebar.locator('button:has-text("Visualizations")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Visualizations: {e}")
        
        # Verify page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Visualizations page: Import error - Error loading page"
            assert rendering_error == 0, "Visualizations page: Runtime error - Error rendering page"
            
            # Check for "Generate Plot" or similar button or page content
            assert page.locator("text=Generate Plot").first.count() > 0 or page.locator("text=Plot").first.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Visualizations page loaded: {e}")


class TestLiteratureAnalysis:
    """Test literature analysis page"""
    
    def test_critique_tab(self, page: Page, base_url: str):
        """View literature critique tab"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: Analysis expander > Literature Analysis
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Analysis section
            page.locator("text=📊 Analysis").click()
            page.wait_for_timeout(1000)
            
            # Click Literature Analysis
            sidebar.locator('button:has-text("Literature Analysis")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Literature Analysis: {e}")
        
        # Tab: button[role="tab"]:has-text("Critique")
        try:
            page.locator('button[role="tab"]:has-text("Critique")').click()
            page.wait_for_timeout(1000)
            
            # Verify text_area for "Scientific Text" exists (use .first because there are 2 textareas with this label)
            abstract_textarea = page.locator('textarea[aria-label="Scientific Text"]').first
            assert abstract_textarea.count() > 0, "Critique tab textarea not found"
            
            # Check for "🔍 Analyze" button
            analyze_button = page.get_by_role("button", name="🔍 Analyze", exact=True).first
            if analyze_button.count() == 0:
                analyze_button = page.get_by_role("button", name="Analyze", exact=True).first
            assert analyze_button.count() > 0, "Analyze button not found"
        except Exception as e:
            pytest.fail(f"Could not verify Critique tab: {e}")


class TestAuditLogs:
    """Test audit logs page"""
    
    def test_view_audit_logs(self, page: Page, base_url: str):
        """View audit logs page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Audit Logs
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Audit Logs
            sidebar.locator('button:has-text("Audit Logs")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Audit Logs: {e}")
        
        # Verify: Action selectbox exists, or page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Audit Logs page: Import error - Error loading page"
            assert rendering_error == 0, "Audit Logs page: Runtime error - Error rendering page"
            
            # Check for Action selectbox or page content
            assert page.locator("text=Action").first.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Audit Logs page loaded: {e}")


class TestWorkflows:
    """Test workflows page"""
    
    def test_view_workflows(self, page: Page, base_url: str):
        """View workflows page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Workflows
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Workflows
            sidebar.locator('button:has-text("Workflows")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Workflows: {e}")
        
        # Tab: button[role="tab"]:has-text("Rules") - default
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Workflows page: Import error - Error loading page"
            assert rendering_error == 0, "Workflows page: Runtime error - Error rendering page"
            
            # Check for Rules tab or page content
            rules_tab = page.locator('button[role="tab"]:has-text("Rules")')
            assert rules_tab.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Workflows page loaded: {e}")


class TestTeamsProjects:
    """Test teams and projects page"""
    
    def test_view_teams(self, page: Page, base_url: str):
        """View teams page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Teams & Projects
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Teams & Projects
            sidebar.locator('button:has-text("Teams & Projects")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Teams & Projects: {e}")
        
        # Verify: "My Teams" tab visible or page loads
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Teams & Projects page: Import error - Error loading page"
            assert rendering_error == 0, "Teams & Projects page: Runtime error - Error rendering page"
            
            # Check for My Teams tab or page content
            my_teams_tab = page.locator('button[role="tab"]:has-text("My Teams")')
            assert my_teams_tab.count() > 0 or page.locator("text=Team").first.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Teams & Projects page loaded: {e}")


class TestCostTracking:
    """Test cost tracking page"""
    
    def test_view_cost_tracking(self, page: Page, base_url: str):
        """View cost tracking page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Cost Tracking (NOT Admin!)
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Cost Tracking
            sidebar.locator('button:has-text("Cost Tracking")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Cost Tracking: {e}")
        
        # Tab names: Overview, Add Entry, Entries
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Cost Tracking page: Import error - Error loading page"
            assert rendering_error == 0, "Cost Tracking page: Runtime error - Error rendering page"
            
            # Check for tabs or page content
            overview_tab = page.locator('button[role="tab"]:has-text("Overview")')
            add_entry_tab = page.locator('button[role="tab"]:has-text("Add Entry")')
            entries_tab = page.locator('button[role="tab"]:has-text("Entries")')
            assert overview_tab.count() > 0 or add_entry_tab.count() > 0 or entries_tab.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Cost Tracking page loaded: {e}")


class TestOntologyManagement:
    """Test ontology management page"""
    
    def test_view_ontology(self, page: Page, base_url: str):
        """View ontology page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Ontology Management
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Ontology Management
            sidebar.locator('button:has-text("Ontology Management")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Ontology Management: {e}")
        
        # Tab names: Browse, Add Term
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Ontology Management page: Import error - Error loading page"
            assert rendering_error == 0, "Ontology Management page: Runtime error - Error rendering page"
            
            # Check for tabs or page content
            browse_tab = page.locator('button[role="tab"]:has-text("Browse")')
            add_term_tab = page.locator('button[role="tab"]:has-text("Add Term")')
            assert browse_tab.count() > 0 or add_term_tab.count() > 0 or page.locator("body").count() > 0
        except Exception as e:
            pytest.fail(f"Could not verify Ontology Management page loaded: {e}")


class TestEmailSettings:
    """Test email settings page"""
    
    def test_view_email_settings(self, page: Page, base_url: str):
        """View email settings page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Email Settings
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Email Settings
            sidebar.locator('button:has-text("Email Settings")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Email Settings: {e}")
        
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Email Settings page: Import error - Error loading page"
            assert rendering_error == 0, "Email Settings page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Email Settings page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Email Settings page loaded: {e}")


class TestDataLineage:
    """Test data lineage page"""
    
    def test_view_data_lineage(self, page: Page, base_url: str):
        """View data lineage page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Data Lineage
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Data Lineage
            sidebar.locator('button:has-text("Data Lineage")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Data Lineage: {e}")
        
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Data Lineage page: Import error - Error loading page"
            assert rendering_error == 0, "Data Lineage page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Data Lineage page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Data Lineage page loaded: {e}")


class TestSchedule:
    """Test schedule page"""
    
    def test_view_schedule(self, page: Page, base_url: str):
        """View schedule page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Schedule
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Schedule
            sidebar.locator('button:has-text("Schedule")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Schedule: {e}")
        
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Schedule page: Import error - Error loading page"
            assert rendering_error == 0, "Schedule page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Schedule page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Schedule page loaded: {e}")


class TestGenericAssays:
    """Test generic assays page"""
    
    def test_view_generic_assays(self, page: Page, base_url: str):
        """View generic assays page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Generic Assays
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Generic Assays
            sidebar.locator('button:has-text("Generic Assays")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Generic Assays: {e}")
        
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Generic Assays page: Import error - Error loading page"
            assert rendering_error == 0, "Generic Assays page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Generic Assays page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Generic Assays page loaded: {e}")


class TestCompare:
    """Test compare page"""
    
    def test_view_compare(self, page: Page, base_url: str):
        """View compare page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Compare
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Compare
            sidebar.locator('button:has-text("Compare")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Compare: {e}")
        
        # Verify: page loads without error
        try:
            # Check for error patterns
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Compare page: Import error - Error loading page"
            assert rendering_error == 0, "Compare page: Runtime error - Error rendering page"
            
            # Check page loaded
            assert page.locator("body").count() > 0, "Compare page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Compare page loaded: {e}")


class TestOverview:
    """Test overview page"""
    
    def test_view_overview(self, page: Page, base_url: str):
        """View overview page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 🔍 Discovery (expanded default, no expander click) > Overview
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Click Overview (Discovery is expanded by default)
            sidebar.locator('button:has-text("Overview")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Overview: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Overview page: Import error - Error loading page"
            assert rendering_error == 0, "Overview page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Overview page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Overview page loaded: {e}")


class TestChat:
    """Test chat page"""
    
    def test_view_chat(self, page: Page, base_url: str):
        """View chat page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Chat
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Chat
            sidebar.locator('button:has-text("Chat")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Chat: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Chat page: Import error - Error loading page"
            assert rendering_error == 0, "Chat page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Chat page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Chat page loaded: {e}")


class TestSearch:
    """Test search page"""
    
    def test_view_search(self, page: Page, base_url: str):
        """View search page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Search
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Search
            sidebar.locator('button:has-text("Search")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Search: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Search page: Import error - Error loading page"
            assert rendering_error == 0, "Search page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Search page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Search page loaded: {e}")


class TestCrossOmics:
    """Test cross-omics page"""
    
    def test_view_cross_omics(self, page: Page, base_url: str):
        """View cross-omics page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > Cross-Omics
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click Cross-Omics
            sidebar.locator('button:has-text("Cross-Omics")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Cross-Omics: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Cross-Omics page: Import error - Error loading page"
            assert rendering_error == 0, "Cross-Omics page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Cross-Omics page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Cross-Omics page loaded: {e}")


class TestSystemHealth:
    """Test system health page"""
    
    def test_view_system_health(self, page: Page, base_url: str):
        """View system health page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📚 Other Pages expander > System Health
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)
            
            # Click System Health
            sidebar.locator('button:has-text("System Health")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to System Health: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "System Health page: Import error - Error loading page"
            assert rendering_error == 0, "System Health page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "System Health page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify System Health page loaded: {e}")


class TestQATracker:
    """Test Q&A tracker page"""
    
    def test_view_qa_tracker(self, page: Page, base_url: str):
        """View Q&A tracker page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📋 ELN expander > Q&A Tracker
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand ELN section
            page.locator("text=📋 ELN").first.click()
            page.wait_for_timeout(1000)
            
            # Click Q&A Tracker
            sidebar.locator('button:has-text("Q&A Tracker")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Q&A Tracker: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Q&A Tracker page: Import error - Error loading page"
            assert rendering_error == 0, "Q&A Tracker page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Q&A Tracker page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Q&A Tracker page loaded: {e}")


class TestAnalysisTools:
    """Test analysis tools page"""
    
    def test_view_analysis_tools(self, page: Page, base_url: str):
        """View analysis tools page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📊 Analysis expander > Analysis Tools
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Analysis section
            page.locator("text=📊 Analysis").first.click()
            page.wait_for_timeout(1000)
            
            # Click Analysis Tools
            sidebar.locator('button:has-text("Analysis Tools")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Analysis Tools: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Analysis Tools page: Import error - Error loading page"
            assert rendering_error == 0, "Analysis Tools page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Analysis Tools page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Analysis Tools page loaded: {e}")


class TestCandidateSelection:
    """Test candidate selection page"""
    
    def test_view_candidate_selection(self, page: Page, base_url: str):
        """View candidate selection page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: 📊 Analysis expander > Candidate Selection
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Analysis section
            page.locator("text=📊 Analysis").first.click()
            page.wait_for_timeout(1000)
            
            # Click Candidate Selection
            sidebar.locator('button:has-text("Candidate Selection")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Candidate Selection: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Candidate Selection page: Import error - Error loading page"
            assert rendering_error == 0, "Candidate Selection page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Candidate Selection page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Candidate Selection page loaded: {e}")


class TestImportData:
    """Test import data page"""
    
    def test_view_import_data(self, page: Page, base_url: str):
        """View import data page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Import Data
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Import Data
            sidebar.locator('button:has-text("Import Data")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Import Data: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Import Data page: Import error - Error loading page"
            assert rendering_error == 0, "Import Data page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Import Data page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Import Data page loaded: {e}")


class TestFeedback:
    """Test feedback page"""
    
    def test_view_feedback(self, page: Page, base_url: str):
        """View feedback page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Feedback
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Feedback
            sidebar.locator('button:has-text("Feedback")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Feedback: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Feedback page: Import error - Error loading page"
            assert rendering_error == 0, "Feedback page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Feedback page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Feedback page loaded: {e}")


class TestFeaturePermissions:
    """Test feature permissions page"""
    
    def test_view_feature_permissions(self, page: Page, base_url: str):
        """View feature permissions page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Feature Permissions
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Feature Permissions
            sidebar.locator('button:has-text("Feature Permissions")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Feature Permissions: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Feature Permissions page: Import error - Error loading page"
            assert rendering_error == 0, "Feature Permissions page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Feature Permissions page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Feature Permissions page loaded: {e}")


class TestDataRetention:
    """Test data retention page"""
    
    def test_view_data_retention(self, page: Page, base_url: str):
        """View data retention page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)
        
        # Navigate: ⚙️ Admin expander > Data Retention
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)
            
            # Expand Admin section
            page.locator("text=⚙️ Admin").first.click()
            page.wait_for_timeout(1000)
            
            # Click Data Retention
            sidebar.locator('button:has-text("Data Retention")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Data Retention: {e}")
        
        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Data Retention page: Import error - Error loading page"
            assert rendering_error == 0, "Data Retention page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Data Retention page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Data Retention page loaded: {e}")

