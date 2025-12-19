"""Help chat assistant for the dashboard."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger
from scripts.dashboard.help_content import HELP_TOPICS

logger = get_logger(__name__)


def get_help_response(question: str) -> str:
    """
    Generate a helpful response to a user question about the platform.

    Args:
        question: User's question

    Returns:
        Helpful response string
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    # Build context from help topics
    help_context = []
    for page_name, help_data in HELP_TOPICS.items():
        help_context.append(
            f"**{help_data['title']} ({page_name}):**\n"
            f"Description: {help_data['description']}\n"
            f"Tips: {'; '.join(help_data.get('tips', []))}"
        )

    context_text = "\n\n".join(help_context)

    prompt = (
        "You are a helpful assistant for the Amprenta Multi-Omics Platform. "
        "Answer the user's question using the provided help documentation.\n\n"
        f"Help Documentation:\n{context_text}\n\n"
        f"User Question: {question}\n\n"
        "Provide a clear, concise answer. If the question is about a specific page, "
        "reference the page name and provide relevant tips. If you don't know the answer, "
        "suggest checking the help documentation or contacting support."
    )

    try:
        logger.info("[HELP_CHAT] Generating response for question: %s", question[:50])
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
        )
        response = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        logger.debug("[HELP_CHAT] Generated response: %s", response[:100])
        return response
    except Exception as e:
        logger.error("[HELP_CHAT] Error generating response: %r", e)
        return f"I apologize, but I encountered an error while generating a response. Please try again or check the help documentation directly."


def render_help_chat():
    """Render the help chat interface."""
    st.subheader("💬 Help Chat")
    st.markdown("Ask questions about the platform and get instant help.")

    # Initialize chat history
    if "chat_history" not in st.session_state:
        st.session_state["chat_history"] = []

    # Display chat history
    if st.session_state["chat_history"]:
        st.markdown("### Chat History")
        for i, message in enumerate(st.session_state["chat_history"]):
            role = message["role"]
            content = message["content"]

            if role == "user":
                with st.chat_message("user"):
                    st.write(content)
            else:
                with st.chat_message("assistant"):
                    st.write(content)

    # Chat input
    question = st.chat_input("Ask a question about the platform...")

    if question:
        # Add user message to history
        st.session_state["chat_history"].append({"role": "user", "content": question})

        # Get response
        with st.spinner("Thinking..."):
            response = get_help_response(question)
            st.session_state["chat_history"].append({"role": "assistant", "content": response})

        st.rerun()

    # Clear button
    if st.session_state["chat_history"]:
        if st.button("Clear Chat History", type="secondary"):
            st.session_state["chat_history"] = []
            st.rerun()
