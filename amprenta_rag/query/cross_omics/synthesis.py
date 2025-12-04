"""
LLM synthesis functions for cross-omics summaries.

Handles the actual OpenAI API calls to synthesize cross-omics summaries
from retrieved context chunks.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.prompt_templates import (
    MAX_CONTEXT_CHUNKS,
    get_cross_omics_system_prompt,
    prepare_context_chunks,
)

logger = get_logger(__name__)


def synthesize_cross_omics_summary(
    prompt: str,
    context_chunks: List[str],
    max_chunks: int = MAX_CONTEXT_CHUNKS,
) -> str:
    """
    Use OpenAI to synthesize a cross-omics summary.
    
    Args:
        prompt: User prompt describing what to summarize
        context_chunks: List of context chunk text strings
        max_chunks: Maximum chunks to include in context
        
    Returns:
        Synthesized summary text
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    # Prepare context
    context = prepare_context_chunks(context_chunks, max_chunks)
    
    system_prompt = get_cross_omics_system_prompt()
    user_content = f"{prompt}\n\nRelevant context chunks:\n{context}"
    
    try:
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            temperature=0.3,
        )
        return resp.choices[0].message.content.strip()  # type: ignore[union-attr]
    except Exception as e:
        logger.error(
            "[RAG][CROSS-OMICS] OpenAI API error synthesizing summary: %r",
            e,
        )
        raise

