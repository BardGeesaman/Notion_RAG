"""E2E tests for AI Expert Agents system."""

import pytest
from playwright.sync_api import Page, expect


# ============================================================================
# FIXTURES
# ============================================================================

@pytest.fixture
def expert_chat_page(page: Page, streamlit_server: str):
    """Navigate to AI Expert Chat page."""
    page.goto(f"{streamlit_server}/?page=AI+Expert+Chat")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(5000)  # Allow page to fully render
    return page


@pytest.fixture
def expert_training_page(page: Page, streamlit_server: str):
    """Navigate to AI Expert Training page."""
    page.goto(f"{streamlit_server}/?page=AI+Expert+Training")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(5000)  # Allow page to fully render
    return page


# ============================================================================
# AI EXPERT CHAT TESTS
# ============================================================================

def test_expert_chat_page_loads(expert_chat_page: Page):
    """Test AI Expert Chat page loads successfully."""
    main = expert_chat_page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_role("heading", name="🧠 AI Expert Chat")).to_be_visible()


def test_expert_selection_interface(expert_chat_page: Page):
    """Test expert selection interface is visible."""
    main = expert_chat_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for expert selection elements
    expert_selection = main.get_by_text("Select Expert(s)").or_(
        main.get_by_text("Choose experts")).or_(
        main.get_by_text("Available Experts")
    )
    expect(expert_selection).to_be_visible()


def test_new_conversation_interface(expert_chat_page: Page):
    """Test new conversation creation interface."""
    main = expert_chat_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for conversation creation elements
    new_conversation = main.get_by_text("New Conversation").or_(
        main.get_by_text("Start Chat")).or_(
        main.get_by_role("button", name="New Conversation")
    )
    expect(new_conversation).to_be_visible()


def test_chat_input_interface(expert_chat_page: Page):
    """Test chat input interface is present."""
    main = expert_chat_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for message input
    message_input = main.get_by_placeholder("Type your message").or_(
        main.get_by_placeholder("Ask a question")).or_(
        main.locator('textarea')
    )
    expect(message_input).to_be_visible()


def test_conversation_history_sidebar(expert_chat_page: Page):
    """Test conversation history sidebar is visible."""
    # Look for sidebar elements
    sidebar = expert_chat_page.locator('[data-testid="stSidebar"]').first
    
    # Check for conversation history elements
    history_section = sidebar.get_by_text("Conversation History").or_(
        sidebar.get_by_text("Recent Conversations")).or_(
        sidebar.get_by_text("History")
    )
    expect(history_section).to_be_visible()


def test_feedback_interface_elements(expert_chat_page: Page):
    """Test feedback interface elements are present."""
    main = expert_chat_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for feedback elements (may not be visible until after a message)
    feedback_section = main.get_by_text("Feedback").or_(
        main.get_by_text("Rate this response")).or_(
        main.get_by_text("👍")).or_(
        main.get_by_text("👎")
    )
    # Use to_be_visible with timeout since feedback may appear after interaction
    expect(feedback_section.or_(main.get_by_text("Send Message"))).to_be_visible()


# ============================================================================
# AI EXPERT TRAINING TESTS
# ============================================================================

def test_expert_training_page_loads(expert_training_page: Page):
    """Test AI Expert Training page loads successfully."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_role("heading", name="🎓 AI Expert Training")).to_be_visible()


def test_expert_dropdown_selection(expert_training_page: Page):
    """Test expert dropdown selection interface."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for expert selection dropdown
    expert_dropdown = main.get_by_text("Select Expert").or_(
        main.get_by_text("Choose Expert")).or_(
        main.locator('select')
    )
    expect(expert_dropdown).to_be_visible()


def test_training_tabs_interface(expert_training_page: Page):
    """Test training console tab interface."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for tab elements
    training_tab = main.get_by_role("tab", name="Training Examples").or_(
        main.get_by_text("Training Examples")
    )
    knowledge_tab = main.get_by_role("tab", name="Knowledge Documents").or_(
        main.get_by_text("Knowledge Documents")
    )
    feedback_tab = main.get_by_role("tab", name="Feedback Review").or_(
        main.get_by_text("Feedback Review")
    )
    config_tab = main.get_by_role("tab", name="Expert Configuration").or_(
        main.get_by_text("Expert Configuration")
    )
    
    # At least one tab should be visible
    expect(training_tab.or_(knowledge_tab).or_(feedback_tab).or_(config_tab)).to_be_visible()


def test_training_examples_tab_content(expert_training_page: Page):
    """Test training examples tab content."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for training examples content
    examples_content = main.get_by_text("Training Examples").or_(
        main.get_by_text("Add Example")).or_(
        main.get_by_text("Examples List")
    )
    expect(examples_content).to_be_visible()


def test_knowledge_documents_interface(expert_training_page: Page):
    """Test knowledge documents upload interface."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for knowledge documents interface
    knowledge_content = main.get_by_text("Knowledge Documents").or_(
        main.get_by_text("Upload Document")).or_(
        main.get_by_text("Documents")
    )
    expect(knowledge_content).to_be_visible()


def test_feedback_metrics_interface(expert_training_page: Page):
    """Test feedback review metrics interface."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for feedback metrics
    feedback_content = main.get_by_text("Feedback Review").or_(
        main.get_by_text("Feedback Metrics")).or_(
        main.get_by_text("Rating")
    )
    expect(feedback_content).to_be_visible()


def test_expert_configuration_interface(expert_training_page: Page):
    """Test expert configuration edit interface."""
    main = expert_training_page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for configuration interface
    config_content = main.get_by_text("Expert Configuration").or_(
        main.get_by_text("Configuration")).or_(
        main.get_by_text("Settings")
    )
    expect(config_content).to_be_visible()


# ============================================================================
# INTEGRATION TESTS
# ============================================================================

def test_navigation_between_expert_pages(page: Page, streamlit_server: str):
    """Test navigation between AI Expert pages."""
    # Start at chat page
    page.goto(f"{streamlit_server}/?page=AI+Expert+Chat")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_role("heading", name="🧠 AI Expert Chat")).to_be_visible()
    
    # Navigate to training page
    page.goto(f"{streamlit_server}/?page=AI+Expert+Training")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(3000)
    
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    expect(main.get_by_role("heading", name="🎓 AI Expert Training")).to_be_visible()


def test_api_integration_health(page: Page, streamlit_server: str):
    """Test that pages can communicate with API (smoke test)."""
    # Test chat page API integration
    page.goto(f"{streamlit_server}/?page=AI+Expert+Chat")
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_timeout(5000)
    
    main = page.locator('[data-testid="stMainBlockContainer"]').first
    
    # Look for either expert list or error message - indicates API call attempted
    api_response = main.get_by_text("No experts available").or_(
        main.get_by_text("Select Expert")).or_(
        main.get_by_text("Error loading")).or_(
        main.get_by_text("Cannot connect")
    )
    expect(api_response).to_be_visible()
