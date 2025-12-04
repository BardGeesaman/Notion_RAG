"""
Answer synthesis for RAG queries.

This module provides functions for synthesizing answers using OpenAI.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def synthesize_answer(user_query: str, chunks: List[str]) -> str:
    """
    Use OpenAI to synthesize an answer given the query + context chunks.

    Args:
        user_query: User's query text
        chunks: List of context chunk texts from Notion

    Returns:
        Synthesized answer string from OpenAI
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    context = "\n\n---\n\n".join(chunks[:8])  # cap context for now

    system_prompt = (
        "You are an assistant helping interpret Amprenta's internal RAG context.\n"
        "Use only the provided context chunks to answer the user's question.\n"
        "If the context is insufficient, say so explicitly."
    )

    user_content = (
        f"User question:\n{user_query}\n\n" f"Relevant context chunks:\n{context}"
    )

    try:
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            temperature=0.2,
        )
    except Exception as e:
        logger.error("[RAG] OpenAI API error synthesizing answer: %r", e)
        raise
    return resp.choices[0].message.content.strip()  # type: ignore[union-attr]

