"""Hypothetical Document Embeddings (HyDE) for improved RAG retrieval."""
from __future__ import annotations

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def generate_hypothetical_answer(query: str) -> str:
    """
    Generate a hypothetical answer to the query for HyDE retrieval.
    
    Args:
        query: User's query text
        
    Returns:
        Hypothetical answer string (2-3 sentences)
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    prompt = f"Write a 2-3 sentence factual answer to this question: {query}"
    
    try:
        logger.info("[HyDE] Generating hypothetical answer for query: %s", query[:50])
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
        )
        answer = resp.choices[0].message.content.strip()  # type: ignore[union-attr]
        logger.debug("[HyDE] Generated answer: %s", answer[:100])
        return answer
    except Exception as e:
        logger.error("[HyDE] OpenAI API error generating hypothetical answer: %r", e)
        raise
