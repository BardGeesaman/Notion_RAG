"""AI Expert Chat Dashboard - Interactive chat with domain experts."""

from __future__ import annotations

import streamlit as st
import requests
from datetime import datetime
from uuid import UUID
import time

st.set_page_config(page_title="AI Expert Chat", page_icon="🧠", layout="wide")

API_BASE = "http://localhost:8000/api/v1/experts"

# Initialize session state
if "current_conversation_id" not in st.session_state:
    st.session_state.current_conversation_id = None
if "selected_experts" not in st.session_state:
    st.session_state.selected_experts = []
if "chat_messages" not in st.session_state:
    st.session_state.chat_messages = []
if "entity_context" not in st.session_state:
    st.session_state.entity_context = None


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_experts():
    """Fetch available experts."""
    try:
        resp = requests.get(f"{API_BASE}/")
        return resp.json() if resp.ok else []
    except:
        return []


def get_user_conversations():
    """Fetch user's conversation history."""
    try:
        resp = requests.get(f"{API_BASE}/conversations")
        return resp.json() if resp.ok else []
    except:
        return []


def create_new_conversation(expert_ids, title=None):
    """Create new conversation with selected experts."""
    try:
        payload = {
            "expert_ids": expert_ids,
            "title": title or f"Chat with {len(expert_ids)} expert(s)",
        }
        if st.session_state.entity_context:
            payload["entity_context"] = st.session_state.entity_context
        
        resp = requests.post(f"{API_BASE}/conversations", json=payload)
        if resp.ok:
            return resp.json()
    except Exception as e:
        st.error(f"Failed to create conversation: {e}")
    return None


def send_message(conversation_id, message):
    """Send message and get expert response."""
    try:
        payload = {"content": message}
        if st.session_state.entity_context:
            payload["entity_context"] = st.session_state.entity_context
        
        resp = requests.post(f"{API_BASE}/conversations/{conversation_id}/messages", json=payload)
        if resp.ok:
            return resp.json()
    except Exception as e:
        st.error(f"Failed to send message: {e}")
    return None


def get_panel_response(expert_ids, question):
    """Get panel discussion response."""
    try:
        payload = {
            "expert_ids": expert_ids,
            "question": question,
        }
        if st.session_state.entity_context:
            payload["entity_context"] = st.session_state.entity_context
        
        resp = requests.post(f"{API_BASE}/panel", json=payload)
        if resp.ok:
            return resp.json()
    except Exception as e:
        st.error(f"Panel discussion failed: {e}")
    return None


def submit_feedback(message_id, rating, correction=None):
    """Submit feedback on expert response."""
    try:
        payload = {"rating": rating}
        if correction:
            payload["correction"] = correction
        
        resp = requests.post(f"{API_BASE}/messages/{message_id}/feedback", json=payload)
        return resp.ok
    except:
        return False


# ============================================================================
# MAIN LAYOUT
# ============================================================================

st.title("🧠 AI Expert Chat")
st.markdown("Chat with domain experts: medicinal chemists, computational biologists, clinicians, and biostatisticians.")

# Sidebar for expert selection and conversation history
with st.sidebar:
    st.header("🎯 Select Experts")
    
    experts = get_experts()
    expert_options = {}
    
    if experts:
        for expert in experts:
            # Create expert display with emoji and specializations
            emoji = expert.get("avatar_emoji", "👨‍🔬")
            display_name = f"{emoji} {expert['name']}"
            specializations = ", ".join(expert.get("specializations", [])[:2])
            
            is_selected = st.checkbox(
                display_name,
                key=f"expert_{expert['id']}",
                help=f"{expert['role']} - {specializations}"
            )
            
            if is_selected:
                expert_options[expert['id']] = expert
        
        st.session_state.selected_experts = list(expert_options.keys())
    else:
        st.warning("No experts available. Check API connection.")
    
    # Panel mode indicator
    if len(st.session_state.selected_experts) > 1:
        st.info(f"🎭 Panel Mode: {len(st.session_state.selected_experts)} experts")
    elif len(st.session_state.selected_experts) == 1:
        selected_expert = expert_options[st.session_state.selected_experts[0]]
        st.success(f"💬 Chat with {selected_expert['name']}")
    else:
        st.info("Select expert(s) to start chatting")
    
    st.markdown("---")
    
    # Entity context selector
    st.header("🔗 Entity Context")
    entity_type = st.selectbox(
        "Link to entity (optional)",
        ["None", "Compound", "Dataset", "Experiment"],
        index=0
    )
    
    if entity_type != "None":
        entity_id = st.text_input("Entity ID (UUID)", placeholder="e.g., 123e4567-...")
        if entity_id:
            try:
                st.session_state.entity_context = {
                    "entity_type": entity_type,
                    "entity_id": entity_id
                }
                st.success(f"Context: {entity_type}")
            except:
                st.error("Invalid UUID format")
        else:
            st.session_state.entity_context = None
    else:
        st.session_state.entity_context = None
    
    st.markdown("---")
    
    # Conversation history
    st.header("📚 Recent Conversations")
    conversations = get_user_conversations()
    
    if conversations:
        for conv in conversations[:10]:  # Show last 10
            expert_names = ", ".join(conv.get("expert_names", []))
            title = conv["title"][:30] + "..." if len(conv["title"]) > 30 else conv["title"]
            
            if st.button(f"💬 {title}", key=f"conv_{conv['id']}", help=f"Experts: {expert_names}"):
                st.session_state.current_conversation_id = conv["id"]
                st.rerun()
    else:
        st.info("No previous conversations")
    
    # New conversation button
    if st.button("➕ New Conversation", type="primary", disabled=not st.session_state.selected_experts):
        if st.session_state.selected_experts:
            conv = create_new_conversation(st.session_state.selected_experts)
            if conv:
                st.session_state.current_conversation_id = conv["id"]
                st.session_state.chat_messages = []
                st.success("New conversation started!")
                st.rerun()


# ============================================================================
# MAIN CHAT AREA
# ============================================================================

# Load current conversation if selected
if st.session_state.current_conversation_id:
    try:
        resp = requests.get(f"{API_BASE}/conversations/{st.session_state.current_conversation_id}")
        if resp.ok:
            conversation = resp.json()
            st.session_state.chat_messages = conversation.get("messages", [])
            
            # Display conversation header
            expert_names = [expert["name"] for expert in conversation.get("experts", [])]
            st.subheader(f"💬 {conversation['title']}")
            st.caption(f"Experts: {', '.join(expert_names)}")
            
            if conversation.get("context_entity_type"):
                st.info(f"Context: {conversation['context_entity_type']} - {conversation.get('context_entity_id', '')}")
        else:
            st.error("Failed to load conversation")
            st.session_state.current_conversation_id = None
    except Exception as e:
        st.error(f"Error loading conversation: {e}")
        st.session_state.current_conversation_id = None

# Chat interface
if st.session_state.selected_experts:
    # Display chat messages
    chat_container = st.container()
    
    with chat_container:
        if st.session_state.chat_messages:
            for msg in st.session_state.chat_messages:
                role = msg["role"]
                content = msg["content"]
                expert_name = msg.get("expert_name", "Expert")
                
                if role == "user":
                    # User message
                    with st.chat_message("user", avatar="👤"):
                        st.write(content)
                
                elif role == "assistant":
                    # Expert response
                    avatar = "🧪" if "chemist" in expert_name.lower() else "🧠"
                    with st.chat_message("assistant", avatar=avatar):
                        st.write(f"**{expert_name}:**")
                        st.write(content)
                        
                        # Reasoning if available
                        if msg.get("reasoning"):
                            with st.expander("🤔 Reasoning"):
                                st.write(msg["reasoning"])
                        
                        # Token count
                        if msg.get("token_count"):
                            st.caption(f"Tokens: {msg['token_count']}")
                        
                        # Feedback interface
                        msg_id = msg["id"]
                        col1, col2, col3 = st.columns([1, 2, 1])
                        
                        with col1:
                            if st.button("👍", key=f"up_{msg_id}", help="Helpful"):
                                if submit_feedback(msg_id, 5):
                                    st.success("Feedback sent!")
                        
                        with col2:
                            rating = st.select_slider(
                                "Rate response",
                                options=[1, 2, 3, 4, 5],
                                value=3,
                                format_func=lambda x: "⭐" * x,
                                key=f"rating_{msg_id}"
                            )
                            if st.button("Submit Rating", key=f"submit_{msg_id}"):
                                if submit_feedback(msg_id, rating):
                                    st.success(f"Rated {rating} stars!")
                        
                        with col3:
                            if st.button("👎", key=f"down_{msg_id}", help="Not helpful"):
                                if submit_feedback(msg_id, 1):
                                    st.success("Feedback sent!")
        else:
            st.info("Start a conversation by typing a message below")
    
    # Message input
    st.markdown("---")
    
    # Panel vs Single expert mode
    if len(st.session_state.selected_experts) > 1:
        # Panel mode
        col1, col2 = st.columns([4, 1])
        with col1:
            user_input = st.text_area(
                "Ask the expert panel:",
                placeholder="e.g., What's the best approach for optimizing this compound's ADMET properties?",
                height=100,
                key="panel_input"
            )
        with col2:
            st.markdown("<br>", unsafe_allow_html=True)
            if st.button("🎭 Ask Panel", type="primary", disabled=not user_input):
                if user_input:
                    with st.spinner("Getting panel response..."):
                        panel_result = get_panel_response(st.session_state.selected_experts, user_input)
                        
                        if panel_result:
                            # Display panel summary
                            st.success(f"Panel Response (Confidence: {panel_result['confidence_score']:.1%})")
                            
                            # Show primary recommendation
                            with st.chat_message("assistant", avatar="🎭"):
                                st.write("**Panel Consensus:**")
                                st.write(panel_result["primary_recommendation"])
                            
                            # Show individual expert responses
                            with st.expander("Individual Expert Responses"):
                                for resp in panel_result["expert_responses"]:
                                    st.markdown(f"**{resp['expert_name']} ({resp['expert_role']}):**")
                                    st.write(resp["content"])
                                    st.markdown("---")
                            
                            st.rerun()
    else:
        # Single expert mode
        col1, col2 = st.columns([4, 1])
        with col1:
            user_input = st.text_area(
                "Your message:",
                placeholder="Ask your question...",
                height=100,
                key="single_input"
            )
        with col2:
            st.markdown("<br>", unsafe_allow_html=True)
            if st.button("💬 Send", type="primary", disabled=not user_input):
                if user_input and st.session_state.current_conversation_id:
                    with st.spinner("Getting expert response..."):
                        response = send_message(st.session_state.current_conversation_id, user_input)
                        
                        if response:
                            st.success("Response received!")
                            st.rerun()
                        else:
                            st.error("Failed to get response")
                elif user_input and not st.session_state.current_conversation_id:
                    # Create new conversation
                    conv = create_new_conversation(st.session_state.selected_experts, "New Chat")
                    if conv:
                        st.session_state.current_conversation_id = conv["id"]
                        with st.spinner("Getting expert response..."):
                            response = send_message(st.session_state.current_conversation_id, user_input)
                            if response:
                                st.success("Conversation started!")
                                st.rerun()

else:
    # No experts selected
    st.info("👈 Select one or more experts from the sidebar to start chatting")
    
    # Show expert showcase
    experts = get_experts()
    if experts:
        st.markdown("### Available Experts")
        
        cols = st.columns(2)
        for i, expert in enumerate(experts):
            with cols[i % 2]:
                with st.container(border=True):
                    emoji = expert.get("avatar_emoji", "👨‍🔬")
                    st.markdown(f"## {emoji} {expert['name']}")
                    st.markdown(f"**{expert['role']}**")
                    
                    specializations = expert.get("specializations", [])
                    if specializations:
                        st.write("**Specializations:**")
                        for spec in specializations[:4]:  # Show first 4
                            st.write(f"• {spec}")
                    
                    st.caption(f"Version: {expert.get('prompt_version', '1.0')}")


# ============================================================================
# CONVERSATION MANAGEMENT
# ============================================================================

# Conversation controls (if active conversation)
if st.session_state.current_conversation_id:
    st.markdown("---")
    col1, col2, col3 = st.columns([2, 2, 1])
    
    with col1:
        if st.button("🔄 Refresh Messages"):
            # Reload conversation
            try:
                resp = requests.get(f"{API_BASE}/conversations/{st.session_state.current_conversation_id}")
                if resp.ok:
                    conversation = resp.json()
                    st.session_state.chat_messages = conversation.get("messages", [])
                    st.rerun()
            except Exception as e:
                st.error(f"Failed to refresh: {e}")
    
    with col2:
        if st.button("📊 Expert Stats"):
            # Show expert statistics
            if st.session_state.selected_experts:
                expert_id = st.session_state.selected_experts[0]
                try:
                    resp = requests.get(f"{API_BASE}/{expert_id}/stats")
                    if resp.ok:
                        stats = resp.json()
                        with st.expander("Expert Performance"):
                            st.metric("Messages", stats.get("message_count", 0))
                            st.metric("Avg Rating", f"{stats.get('average_rating', 0):.1f}/5.0" if stats.get('average_rating') else "N/A")
                            st.metric("Knowledge Docs", stats.get("knowledge_doc_count", 0))
                except:
                    st.error("Failed to load stats")
    
    with col3:
        if st.button("🗑️ End Chat"):
            if st.session_state.current_conversation_id:
                try:
                    resp = requests.delete(f"{API_BASE}/conversations/{st.session_state.current_conversation_id}")
                    if resp.ok:
                        st.session_state.current_conversation_id = None
                        st.session_state.chat_messages = []
                        st.success("Conversation ended")
                        st.rerun()
                except:
                    st.error("Failed to end conversation")


# ============================================================================
# TIPS AND HELP
# ============================================================================

with st.expander("💡 Tips for Better Conversations"):
    st.markdown("""
    **Getting the best responses:**
    - Be specific about your question or problem
    - Provide context about your research area
    - Link to specific compounds/datasets for targeted advice
    - Use panel mode for complex questions requiring multiple perspectives
    
    **Expert Specializations:**
    - **Med Chemist**: SAR analysis, ADMET optimization, lead design
    - **Comp Biologist**: Multi-omics, pathway analysis, target validation  
    - **Clinician**: Trial design, biomarkers, translational medicine
    - **Biostatistician**: Power analysis, dose-response, experimental design
    
    **Panel Discussions:**
    - Select 2+ experts for collaborative responses
    - Great for complex questions with multiple angles
    - Experts can build on each other's insights
    - Consensus scoring shows agreement level
    """)


# ============================================================================
# REGISTER IN PAGE_REGISTRY
# ============================================================================

def main():
    """Entry point for page registry."""
    return  # Page renders on import
