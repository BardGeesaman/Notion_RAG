"""
Answer synthesis for RAG queries.

This module provides functions for synthesizing answers using OpenAI.
"""

from __future__ import annotations

from typing import List, Tuple, Dict, Any

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.rag.models import Citation
from amprenta_rag.llm.model_registry import get_model_client, AVAILABLE_MODELS

logger = get_logger(__name__)


def synthesize_answer(user_query: str, chunks: List[str], model: str = "gpt-4o") -> str:
    """
    Use OpenAI to synthesize an answer given the query + context chunks.

    Args:
        user_query: User's query text
        chunks: List of context chunk texts from Notion
        model: Model name to use (default: "gpt-4o")

    Returns:
        Synthesized answer string from OpenAI
    """
    client = get_model_client(model)
    actual_model_name = AVAILABLE_MODELS[model]["name"]

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
            model=actual_model_name,
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


def build_citation_context(chunks: List[str], metadata_list: List[Dict[str, Any]]) -> Tuple[str, List[Citation]]:
    """Build numbered context and citation objects."""
    citations: List[Citation] = []
    numbered_chunks: List[str] = []
    for i, (chunk, meta) in enumerate(zip(chunks, metadata_list), 1):
        numbered_chunks.append(f"[{i}] {chunk}")
        citations.append(
            Citation(
                number=i,
                chunk_id=meta.get("chunk_id", meta.get("id", f"chunk_{i}")),
                source_type=meta.get("source_type"),
                title=meta.get("title"),
                dataset_name=meta.get("dataset_name"),
                experiment_name=meta.get("experiment_name"),
                url=meta.get("url"),
            )
        )
    context = "\n\n".join(numbered_chunks)
    return context, citations


def synthesize_answer_with_citations(
    user_query: str,
    chunks: List[str],
    metadata_list: List[Dict[str, Any]],
    model: str = "gpt-4o",
) -> Tuple[str, List[Citation]]:
    """Generate answer with inline citations."""
    context, citations = build_citation_context(chunks, metadata_list)

    prompt = (
        "Answer the following question using ONLY the provided sources.\n"
        "Include inline citations [1], [2], etc. when referencing information from sources.\n\n"
        f"Sources:\n{context}\n\n"
        f"Question: {user_query}\n\n"
        "Answer (cite sources with [1], [2], etc.):"
    )

    client = get_model_client(model)
    actual_model_name = AVAILABLE_MODELS[model]["name"]
    response = client.chat.completions.create(
        model=actual_model_name,
        messages=[{"role": "user", "content": prompt}],
        temperature=0.2,
    )
    answer = response.choices[0].message.content
    return answer, citations

