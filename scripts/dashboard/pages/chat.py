from datetime import datetime
from uuid import uuid4

import streamlit as st

from amprenta_rag.agent.chat_agent import run_chat_turn
from amprenta_rag.agent.chat_types import ChatSessionState
from scripts.dashboard.auth import require_auth


def get_session() -> ChatSessionState:
    if "chat_session" not in st.session_state:
        st.session_state.chat_session = ChatSessionState(id=uuid4(), created_at=datetime.utcnow(), turns=[])
    return st.session_state.chat_session


def render_chat_page():
    require_auth()
    st.title("💬 Amprenta AI Assistant")

    # Info box with instructions
    if "show_instructions" not in st.session_state:
        st.session_state.show_instructions = True

    if st.session_state.show_instructions:
        with st.expander("ℹ️ How to use the chat", expanded=False):
            st.markdown(
                """
            **Available commands:**
            - **Freeform questions:** "What do we know about ceramide in ALS?"
            - **Dataset summary:** "Tell me about dataset: <UUID>"
            - **Program summary:** "Program summary for: <program-name>"
            - **Feature info:** "Feature: TP53"
            - **Signature info:** "What is signature: <ID>?"
            - **Help:** Type "help" to see all commands

            💡 Tip: Get entity IDs from the Datasets, Programs, Features, or Signatures pages
            """
            )

    session = get_session()

    # Display chat history in message bubbles
    if len(session.turns) > 0:
        st.markdown("### 💬 Conversation")
        for i, turn in enumerate(session.turns):
            for msg in turn.messages:
                if msg.role == "user":
                    # User message (right-aligned, blue)
                    st.markdown(
                        f"""
                    <div style="background-color: #1f77b4; padding: 10px; border-radius: 10px; margin: 5px 0; color: white;">
                        <strong>You:</strong> {msg.content}
                    </div>
                    """,
                        unsafe_allow_html=True,
                    )
                else:
                    # Assistant message (left-aligned, gray)
                    st.markdown(
                        f"""
                    <div style="background-color: #2d3748; padding: 10px; border-radius: 10px; margin: 5px 0; color: white;">
                        <strong>Assistant:</strong><br>{msg.content}
                    </div>
                    """,
                        unsafe_allow_html=True,
                    )
        st.markdown("---")
    else:
        st.info("👋 Welcome! Ask me anything about your multi-omics data. Type 'help' to see what I can do.")

    # Input area at bottom
    st.markdown("### Ask a question")
    user_input = st.text_area(
        "Your question",
        key="chat_input",
        height=100,
        placeholder="e.g. What do we know about ceramide in ALS?",
        label_visibility="collapsed",
    )

    col1, col2 = st.columns([3, 1])
    with col1:
        if st.button("Send", type="primary", use_container_width=True):
            text = user_input.strip()
            if not text:
                st.warning("Please enter a question.")
            else:
                with st.spinner("🤔 Thinking..."):
                    session, answer = run_chat_turn(session, text)
                    # Update session state
                    st.session_state.chat_session = session
                    # Rerun to refresh the page and clear input
                    st.rerun()
    with col2:
        if st.button("🗑️ Clear", type="secondary", use_container_width=True):
            st.session_state.chat_session = ChatSessionState(id=uuid4(), created_at=datetime.utcnow(), turns=[])
            st.rerun()


if __name__ == "__main__":
    render_chat_page()
