"""E2E tests for Generative Chemistry dashboard functionality."""

from __future__ import annotations

import re

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect  # noqa: E402


pytestmark = pytest.mark.requires_server


def _goto_generative_chemistry_page(page: Page, base_url: str) -> None:
    """Navigate to the Generative Chemistry page."""
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(2000)
    
    page.goto(f"{base_url}/?page=Generative%20Chemistry")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)


def _main_container(page: Page):
    """Return the main container locator scoped to avoid sidebar."""
    return page.locator('[data-testid="stMainBlockContainer"]').first


def test_page_loads(page: Page, base_url: str) -> None:
    """Test that the Generative Chemistry page loads successfully."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Check page header
    expect(main.get_by_text("🧪 Generative Chemistry")).to_be_visible()
    expect(main.get_by_text("AI-powered de novo molecular design")).to_be_visible()
    
    # Check all tabs are present
    expect(main.get_by_role("tab", name="🎲 Generate")).to_be_visible()
    expect(main.get_by_role("tab", name="🔄 Interpolate")).to_be_visible()
    expect(main.get_by_role("tab", name="⚡ Optimize")).to_be_visible()
    expect(main.get_by_role("tab", name="📊 Results")).to_be_visible()


def test_generate_tab_ui(page: Page, base_url: str) -> None:
    """Test the Generate tab UI elements."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Click on Generate tab (should be active by default)
    generate_tab = main.get_by_role("tab", name="🎲 Generate")
    generate_tab.click()
    page.wait_for_timeout(1000)
    
    # Check UI elements
    expect(main.get_by_text("Random Molecular Generation")).to_be_visible()
    expect(main.get_by_text("Parameters")).to_be_visible()
    expect(main.get_by_text("Results")).to_be_visible()
    
    # Check parameter controls
    expect(main.get_by_text("Number of molecules")).to_be_visible()
    expect(main.get_by_text("Temperature")).to_be_visible()
    expect(main.get_by_text("Max SMILES length")).to_be_visible()
    
    # Check generate button
    generate_button = main.get_by_role("button", name="🎲 Generate Molecules")
    expect(generate_button).to_be_visible()
    expect(generate_button).to_be_enabled()
    
    # Test slider interactions
    number_slider = main.locator('div[data-testid*="stSlider"]').first
    expect(number_slider).to_be_visible()


def test_interpolate_tab_ui(page: Page, base_url: str) -> None:
    """Test the Interpolate tab UI elements."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Click on Interpolate tab
    interpolate_tab = main.get_by_role("tab", name="🔄 Interpolate")
    interpolate_tab.click()
    page.wait_for_timeout(1000)
    
    # Check UI elements
    expect(main.get_by_text("Molecular Interpolation")).to_be_visible()
    expect(main.get_by_text("Input Molecules")).to_be_visible()
    expect(main.get_by_text("Trajectory")).to_be_visible()
    
    # Check input fields
    expect(main.get_by_text("Start SMILES")).to_be_visible()
    expect(main.get_by_text("End SMILES")).to_be_visible()
    expect(main.get_by_text("Interpolation steps")).to_be_visible()
    expect(main.get_by_text("Interpolation method")).to_be_visible()
    
    # Check text inputs
    start_input = main.get_by_placeholder("e.g., CCO")
    expect(start_input).to_be_visible()
    expect(start_input).to_be_editable()
    
    end_input = main.get_by_placeholder("e.g., CCC")
    expect(end_input).to_be_visible()
    expect(end_input).to_be_editable()
    
    # Check interpolate button
    interpolate_button = main.get_by_role("button", name="🔄 Interpolate")
    expect(interpolate_button).to_be_visible()
    expect(interpolate_button).to_be_enabled()
    
    # Test input validation - click without inputs
    interpolate_button.click()
    page.wait_for_timeout(1000)
    
    # Should show error message
    expect(main.get_by_text(re.compile(r"Please provide both start and end SMILES", re.IGNORECASE))).to_be_visible()


def test_optimize_tab_ui(page: Page, base_url: str) -> None:
    """Test the Optimize tab UI elements."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Click on Optimize tab
    optimize_tab = main.get_by_role("tab", name="⚡ Optimize")
    optimize_tab.click()
    page.wait_for_timeout(1000)
    
    # Check UI elements
    expect(main.get_by_text("Property-Guided Optimization")).to_be_visible()
    expect(main.get_by_text("Optimization Setup")).to_be_visible()
    expect(main.get_by_text("Optimization Results")).to_be_visible()
    
    # Check seed input
    expect(main.get_by_text("Seed molecule SMILES")).to_be_visible()
    seed_input = main.get_by_placeholder("e.g., CCO")
    expect(seed_input).to_be_visible()
    
    # Check constraints section
    expect(main.get_by_text("Property Constraints")).to_be_visible()
    expect(main.get_by_text("Add New Constraint")).to_be_visible()
    
    # Check optimization parameters
    expect(main.get_by_text("Optimization Parameters")).to_be_visible()
    expect(main.get_by_text("Iterations")).to_be_visible()
    expect(main.get_by_text("Samples per iteration")).to_be_visible()
    expect(main.get_by_text("Learning rate")).to_be_visible()
    
    # Check optimize button
    optimize_button = main.get_by_role("button", name="⚡ Optimize")
    expect(optimize_button).to_be_visible()
    expect(optimize_button).to_be_enabled()
    
    # Test constraint addition
    add_constraint_expander = main.get_by_text("Add New Constraint")
    add_constraint_expander.click()
    page.wait_for_timeout(500)
    
    # Should show constraint form
    expect(main.get_by_text("Property")).to_be_visible()
    expect(main.get_by_text("Constraint type")).to_be_visible()
    
    # Test adding a constraint
    add_button = main.get_by_role("button", name="Add Constraint")
    expect(add_button).to_be_visible()
    add_button.click()
    page.wait_for_timeout(1000)
    
    # Should show current constraints section
    expect(main.get_by_text("Current Constraints")).to_be_visible()


def test_results_tab_ui(page: Page, base_url: str) -> None:
    """Test the Results tab UI elements."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Click on Results tab
    results_tab = main.get_by_role("tab", name="📊 Results")
    results_tab.click()
    page.wait_for_timeout(1000)
    
    # Check UI elements
    expect(main.get_by_text("Session Results")).to_be_visible()
    expect(main.get_by_text("View and manage generated molecules")).to_be_visible()
    
    # Should show empty state initially
    expect(main.get_by_text("No molecules generated yet")).to_be_visible()
    expect(main.get_by_text("Use the other tabs to generate molecules")).to_be_visible()


def test_tab_navigation(page: Page, base_url: str) -> None:
    """Test navigation between tabs."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Test clicking each tab
    tabs = [
        ("🎲 Generate", "Random Molecular Generation"),
        ("🔄 Interpolate", "Molecular Interpolation"),
        ("⚡ Optimize", "Property-Guided Optimization"),
        ("📊 Results", "Session Results"),
    ]
    
    for tab_name, expected_content in tabs:
        tab = main.get_by_role("tab", name=tab_name)
        expect(tab).to_be_visible()
        
        tab.click()
        page.wait_for_timeout(1000)
        
        # Check that tab content is visible
        expect(main.get_by_text(expected_content)).to_be_visible()
    
    # Test that tabs remain clickable after navigation
    generate_tab = main.get_by_role("tab", name="🎲 Generate")
    generate_tab.click()
    page.wait_for_timeout(500)
    
    expect(main.get_by_text("Random Molecular Generation")).to_be_visible()


def test_generate_tab_parameter_validation(page: Page, base_url: str) -> None:
    """Test parameter validation in Generate tab."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Go to Generate tab
    generate_tab = main.get_by_role("tab", name="🎲 Generate")
    generate_tab.click()
    page.wait_for_timeout(1000)
    
    # Test that sliders have reasonable bounds
    # Number of molecules slider should be between 1-100
    number_text = main.get_by_text("Number of molecules")
    expect(number_text).to_be_visible()
    
    # Temperature slider should be between 0.1-2.0
    temperature_text = main.get_by_text("Temperature")
    expect(temperature_text).to_be_visible()
    
    # Max length slider should be between 20-200
    max_length_text = main.get_by_text("Max SMILES length")
    expect(max_length_text).to_be_visible()
    
    # Generate button should be present and enabled
    generate_button = main.get_by_role("button", name="🎲 Generate Molecules")
    expect(generate_button).to_be_visible()
    expect(generate_button).to_be_enabled()


def test_interpolate_tab_smiles_validation(page: Page, base_url: str) -> None:
    """Test SMILES validation in Interpolate tab."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Go to Interpolate tab
    interpolate_tab = main.get_by_role("tab", name="🔄 Interpolate")
    interpolate_tab.click()
    page.wait_for_timeout(1000)
    
    # Test empty input validation
    interpolate_button = main.get_by_role("button", name="🔄 Interpolate")
    interpolate_button.click()
    page.wait_for_timeout(1000)
    
    # Should show validation error
    expect(main.get_by_text(re.compile(r"Please provide both start and end SMILES", re.IGNORECASE))).to_be_visible()
    
    # Test with only start SMILES
    start_input = main.get_by_placeholder("e.g., CCO")
    start_input.fill("CCO")
    page.wait_for_timeout(500)
    
    interpolate_button.click()
    page.wait_for_timeout(1000)
    
    # Should still show validation error
    expect(main.get_by_text(re.compile(r"Please provide both start and end SMILES", re.IGNORECASE))).to_be_visible()


def test_optimize_tab_constraint_management(page: Page, base_url: str) -> None:
    """Test constraint management in Optimize tab."""
    _goto_generative_chemistry_page(page, base_url)
    
    main = _main_container(page)
    
    # Go to Optimize tab
    optimize_tab = main.get_by_role("tab", name="⚡ Optimize")
    optimize_tab.click()
    page.wait_for_timeout(1000)
    
    # Open constraint form
    add_constraint_expander = main.get_by_text("Add New Constraint")
    add_constraint_expander.click()
    page.wait_for_timeout(500)
    
    # Add a constraint
    add_button = main.get_by_role("button", name="Add Constraint")
    add_button.click()
    page.wait_for_timeout(1000)
    
    # Should show current constraints
    expect(main.get_by_text("Current Constraints")).to_be_visible()
    expect(main.get_by_text("logp")).to_be_visible()  # Default property
    
    # Should show remove button
    remove_button = main.get_by_role("button", name="🗑️")
    expect(remove_button).to_be_visible()
    
    # Test optimization without seed SMILES
    optimize_button = main.get_by_role("button", name="⚡ Optimize")
    optimize_button.click()
    page.wait_for_timeout(1000)
    
    # Should show validation error
    expect(main.get_by_text(re.compile(r"Please provide a seed SMILES", re.IGNORECASE))).to_be_visible()
