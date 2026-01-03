"""AI Expert Chat Dashboard - Interactive chat with domain experts."""

from __future__ import annotations

import streamlit as st
import requests
from datetime import datetime
from uuid import UUID
import time

API_BASE = "http://localhost:8000/api/v1/experts"


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


def get_conversation_messages(conversation_id):
    """Fetch messages for a conversation."""
    try:
        resp = requests.get(f"{API_BASE}/conversations/{conversation_id}/messages")
        return resp.json() if resp.ok else []
    except:
        return []


def record_feedback(message_id, rating, correction=None, tags=None):
    """Record user feedback on expert response."""
    try:
        payload = {
            "rating": rating,
            "correction": correction,
            "tags": tags or []
        }
        resp = requests.post(f"{API_BASE}/feedback/{message_id}", json=payload)
        return resp.ok
    except:
        return False


def display_message(msg, show_feedback=True):
    """Display a chat message with optional feedback controls."""
    timestamp = datetime.fromisoformat(msg.get("timestamp", "").replace("Z", "+00:00"))
    
    if msg["role"] == "user":
        with st.chat_message("user"):
            st.write(msg["content"])
            st.caption(f"You • {timestamp.strftime('%H:%M')}")
    else:
        with st.chat_message("assistant", avatar="🧠"):
            st.write(msg["content"])
            
            # Expert info
            expert_name = msg.get("expert_name", "Expert")
            st.caption(f"{expert_name} • {timestamp.strftime('%H:%M')}")
            
            # Inline feedback
            if show_feedback and msg.get("id"):
                with st.expander("💬 Feedback", expanded=False):
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        rating = st.select_slider(
                            "Rating",
                            options=[1, 2, 3, 4, 5],
                            value=3,
                            key=f"rating_{msg['id']}"
                        )
                    
                    with col2:
                        correction = st.text_area(
                            "Correction/Comment",
                            placeholder="Optional feedback...",
                            key=f"correction_{msg['id']}",
                            height=60
                        )
                    
                    tags = st.multiselect(
                        "Tags",
                        options=["accuracy", "relevance", "clarity", "completeness", "conciseness"],
                        key=f"tags_{msg['id']}"
                    )
                    
                    if st.button("Submit Feedback", key=f"feedback_{msg['id']}"):
                        if record_feedback(msg["id"], rating, correction, tags):
                            st.success("Feedback recorded!")
                        else:
                            st.error("Failed to record feedback")


def main():
    """Entry point for page registry."""
    st.set_page_config(page_title="AI Expert Chat", page_icon="🧠", layout="wide")

    # Initialize session state
    if "current_conversation_id" not in st.session_state:
        st.session_state.current_conversation_id = None
    if "selected_experts" not in st.session_state:
        st.session_state.selected_experts = []
    if "chat_messages" not in st.session_state:
        st.session_state.chat_messages = []
    if "entity_context" not in st.session_state:
        st.session_state.entity_context = None

    st.title("🧠 AI Expert Chat")
    st.markdown("Chat with specialized AI experts for domain-specific guidance and insights.")

    # ============================================================================
    # SIDEBAR: EXPERT SELECTION & CONVERSATION HISTORY
    # ============================================================================

    with st.sidebar:
        st.header("🎯 Expert Selection")
        
        # Load available experts
        experts = get_experts()
        
        if experts:
            expert_options = {f"{exp['name']} ({', '.join(exp['specializations'][:2])})": exp['id'] for exp in experts}
            
            selected_expert_names = st.multiselect(
                "Choose Experts",
                options=list(expert_options.keys()),
                default=[],
                help="Select one or more experts for your conversation"
            )
            
            st.session_state.selected_experts = [expert_options[name] for name in selected_expert_names]
            
            # Show expert details
            if st.session_state.selected_experts:
                st.subheader("Selected Experts")
                for expert in experts:
                    if expert['id'] in st.session_state.selected_experts:
                        with st.expander(expert['name']):
                            st.write(f"**Specializations:** {', '.join(expert['specializations'])}")
                            st.write(f"**Description:** {expert.get('description', 'No description')}")
                            if expert.get('expertise_areas'):
                                st.write(f"**Expertise:** {', '.join(expert['expertise_areas'])}")
        else:
            st.warning("No experts available. Please check the API connection.")
        
        st.divider()
        
        # Entity Context
        st.header("🔗 Entity Context")
        
        entity_type = st.selectbox(
            "Entity Type",
            options=["None", "Compound", "Experiment", "Dataset", "Target"],
            index=0
        )
        
        if entity_type != "None":
            entity_id = st.text_input("Entity ID", placeholder="Enter entity ID...")
            
            if st.button("Set Context"):
                if entity_id:
                    st.session_state.entity_context = {
                        "entity_type": entity_type.lower(),
                        "entity_id": entity_id
                    }
                    st.success(f"Context set: {entity_type} {entity_id}")
                else:
                    st.session_state.entity_context = None
        else:
            st.session_state.entity_context = None
        
        if st.session_state.entity_context:
            st.info(f"Context: {st.session_state.entity_context['entity_type']} {st.session_state.entity_context['entity_id']}")
        
        st.divider()
        
        # Conversation History
        st.header("💬 Recent Conversations")
        
        conversations = get_user_conversations()
        
        if conversations:
            for conv in conversations[:10]:  # Show last 10
                conv_title = conv.get('title', f"Conversation {conv['id']}")
                if st.button(conv_title, key=f"conv_{conv['id']}", use_container_width=True):
                    st.session_state.current_conversation_id = conv['id']
                    st.session_state.chat_messages = get_conversation_messages(conv['id'])
                    st.rerun()
        else:
            st.info("No previous conversations")
        
        # New conversation button
        if st.button("🆕 New Conversation", type="primary", use_container_width=True):
            if st.session_state.selected_experts:
                new_conv = create_new_conversation(st.session_state.selected_experts)
                if new_conv:
                    st.session_state.current_conversation_id = new_conv['id']
                    st.session_state.chat_messages = []
                    st.success("New conversation started!")
                    st.rerun()
            else:
                st.warning("Please select at least one expert first")

    # ============================================================================
    # MAIN CHAT INTERFACE
    # ============================================================================

    if not st.session_state.current_conversation_id:
        # Welcome screen
        st.info("👋 Welcome! Select experts from the sidebar and start a new conversation.")
        
        if experts:
            st.subheader("Available Experts")
            
            cols = st.columns(3)
            for i, expert in enumerate(experts[:6]):  # Show first 6
                with cols[i % 3]:
                    with st.container():
                        st.write(f"**{expert['name']}**")
                        st.write(f"*{', '.join(expert['specializations'][:2])}*")
                        st.write(expert.get('description', '')[:100] + "..." if len(expert.get('description', '')) > 100 else expert.get('description', ''))
                        
                        if expert.get('expertise_areas'):
                            st.caption(f"Expertise: {', '.join(expert['expertise_areas'][:3])}")
        
        # Show visible textarea on welcome screen
        st.text_area("", 
                    placeholder="Start a conversation to begin chatting with experts",
                    disabled=False,  # Make it visible
                    key="welcome_chat_area",
                    height=60)
    else:
        # Active conversation
        st.subheader(f"💬 Conversation")
        
        # Display expert participants
        if st.session_state.selected_experts:
            expert_names = []
            for expert in experts:
                if expert['id'] in st.session_state.selected_experts:
                    expert_names.append(expert['name'])
            
            if expert_names:
                st.caption(f"Chatting with: {', '.join(expert_names)}")
        
        # Context indicator
        if st.session_state.entity_context:
            st.info(f"🔗 Context: {st.session_state.entity_context['entity_type']} {st.session_state.entity_context['entity_id']}")
        
        # Chat messages container
        chat_container = st.container()
        
        with chat_container:
            # Display conversation history
            for msg in st.session_state.chat_messages:
                display_message(msg)
        
        # Chat input (active conversation)
        if prompt := st.chat_input("Ask your experts anything..."):
            # Add user message to display
            user_msg = {
                "role": "user",
                "content": prompt,
                "timestamp": datetime.now().isoformat()
            }
            st.session_state.chat_messages.append(user_msg)
            
            # Display user message immediately
            with chat_container:
                display_message(user_msg, show_feedback=False)
            
            # Send message and get response
            with st.spinner("Experts are thinking..."):
                response = send_message(st.session_state.current_conversation_id, prompt)
                
                if response:
                    # Add expert response to display
                    expert_msg = {
                        "id": response.get("id"),
                        "role": "assistant", 
                        "content": response.get("content", ""),
                        "expert_name": response.get("expert_name", "Expert"),
                        "timestamp": response.get("timestamp", datetime.now().isoformat())
                    }
                    st.session_state.chat_messages.append(expert_msg)
                    
                    # Display expert response
                    with chat_container:
                        display_message(expert_msg)
                else:
                    st.error("Failed to get response from experts")
            
            # Rerun to update the chat
            st.rerun()

    # ============================================================================
    # FOOTER
    # ============================================================================

    st.divider()
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.session_state.current_conversation_id:
            if st.button("🔄 Refresh Messages"):
                st.session_state.chat_messages = get_conversation_messages(st.session_state.current_conversation_id)
                st.rerun()
    
    with col2:
        if st.session_state.current_conversation_id:
            if st.button("🗑️ Clear Chat"):
                st.session_state.chat_messages = []
                st.rerun()
    
    with col3:
        if st.button("ℹ️ Help"):
            st.info("""
            **How to use AI Expert Chat:**
            1. Select one or more experts from the sidebar
            2. Optionally set entity context for domain-specific questions
            3. Start a new conversation or continue an existing one
            4. Ask questions and get expert responses
            5. Provide feedback to help improve expert responses
            
            **Tips:**
            - Use entity context for questions about specific compounds, experiments, etc.
            - Select multiple experts for different perspectives
            - Use the feedback system to improve future responses
            """)

    # Global chat input (always visible)
    if st.session_state.current_conversation_id:
        placeholder_text = "Ask your experts anything..."
    else:
        placeholder_text = "Start a conversation to begin chatting with experts"
    
    if prompt := st.chat_input(placeholder_text):
        if st.session_state.current_conversation_id:
            # Handle active conversation
            user_msg = {
                "role": "user", 
                "content": prompt,
                "timestamp": datetime.now().isoformat()
            }
            st.session_state.chat_messages.append(user_msg)
            st.rerun()
        else:
            st.warning("Please start a conversation first by selecting experts and clicking 'New Conversation'")
