"""
Main RAG query orchestration.

This module provides the main query functions that orchestrate the complete
RAG pipeline: Pinecone queries, match processing, chunk collection, and answer synthesis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
from amprenta_rag.query.semantic_cache import get_semantic_cache

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.bm25_search import bm25_search, reciprocal_rank_fusion
from amprenta_rag.query.pinecone_query import build_meta_filter, query_pinecone
from amprenta_rag.query.rag.chunk_collection import collect_chunks
from amprenta_rag.rag.hybrid_chunk_collection import collect_hybrid_chunks
from amprenta_rag.query.rag.match_processing import filter_matches, summarize_match
from amprenta_rag.query.rag.models import RAGQueryResult, Citation
from amprenta_rag.query.rag.synthesis import synthesize_answer, synthesize_answer_with_citations
from amprenta_rag.query.reranker import get_reranker
from amprenta_rag.query.hyde import generate_hypothetical_answer
from amprenta_rag.query.hallucination import check_groundedness
from amprenta_rag.query.evaluation import evaluate_rag_response, EvalResult
from amprenta_rag.query.agent import agentic_rag

logger = get_logger(__name__)


def query_rag(
    user_query: str,
    *,
    top_k: int = 10,
    source_types: Optional[List[str]] = None,
    tag: Optional[str] = None,
    generate_answer: bool = True,
    disease: Optional[str] = None,
    target: Optional[str] = None,
    lipid: Optional[str] = None,
    signature: Optional[str] = None,
    use_postgres: bool = True,
    use_hybrid: bool = False,
    hybrid_alpha: float = 0.5,
    use_rerank: bool = False,
    rerank_model: Optional[str] = None,
    use_cache: bool = True,
    use_hyde: bool = False,
    check_hallucination: bool = False,
    evaluate: bool = False,
    use_agent: bool = False,
) -> RAGQueryResult:
    """
    High-level API: run a complete RAG query and get structured results.

    This function orchestrates the full RAG pipeline:
    1. Builds metadata filter from query constraints
    2. Queries Pinecone for matches (with source type filtering)
    3. Optionally runs BM25 / Postgres full‑text search and fuses results
       with Reciprocal Rank Fusion when ``use_hybrid=True``
    4. Summarizes and filters matches
    5. Retrieves chunk text from Postgres (and legacy Notion snippets as fallback)
    6. Synthesizes answer using OpenAI (if enabled)

    Args:
        user_query: User's query text
        top_k: Number of results to retrieve from Pinecone
        source_types: Optional list of source types to filter by (e.g., ["Literature", "Experiment"])
        tag: Optional tag substring filter
        generate_answer: Whether to synthesize an answer (default: True)
        disease: Optional disease filter
        target: Optional molecular target filter
        lipid: Optional lipid filter (canonical ID or raw label)
        signature: Optional lipid signature filter
        use_rerank: Whether to rerank results with a cross-encoder
        rerank_model: Optional cross-encoder model name to use for reranking

    Returns:
        RAGQueryResult with matches, filtered matches, context chunks, and answer
    """
    # Check semantic cache
    if use_cache and generate_answer:
        cache = get_semantic_cache()
        cached = cache.get(user_query)
        if cached is not None:
            logger.info("[RAG] Returning cached result")
            return cached
    
    # Agentic multi-step flow
    if use_agent:
        logger.info("[RAG] Using agentic multi-step flow")
        agent_result = agentic_rag(user_query, max_steps=3)
        return RAGQueryResult(
            query=user_query,
            answer=agent_result.answer,
            matches=[],
            filtered_matches=[],
            context_chunks=[],
            citations=[],
        )
    
    # Generate hypothetical answer for HyDE if enabled
    if use_hyde:
        logger.info("[RAG][HYDE] Generating hypothetical answer")
        search_query = generate_hypothetical_answer(user_query)
    else:
        search_query = user_query
    
    logger.info("[RAG] Querying Pinecone (top_k=%d) for: %s", top_k, user_query)
    meta_filter = build_meta_filter(
        disease=disease,
        target=target,
        lipid=lipid,
        signature=signature,
    )

    # 1) Primary vector search (Pinecone)
    raw_matches = query_pinecone(
        search_query,
        top_k=top_k,
        meta_filter=meta_filter,
        source_types=source_types,
    )

    # 2) Optional hybrid search with Postgres BM25
    if use_hybrid:
        logger.info("[RAG][HYBRID] Running hybrid search (vector + BM25, alpha=%.2f)", hybrid_alpha)
        try:
            bm25_results = bm25_search(
                user_query,
                limit=top_k,
                # bm25_search currently supports a single source_type filter
                source_type=source_types[0] if source_types and len(source_types) == 1 else None,
            )
        except Exception as e:
            logger.warning("[RAG][HYBRID] BM25 search failed, falling back to vector-only: %r", e)
            bm25_results = []

        # Normalize Pinecone matches into simple dicts for fusion
        vector_results: List[Dict[str, Any]] = []
        for m in raw_matches:
            # Pinecone client may return objects or dicts
            meta = getattr(m, "metadata", None) or getattr(m, "get", lambda *_: {})(  # type: ignore[attr-defined]
                "metadata", {}
            )
            # Fallback if getattr(get) above returned a callable default
            if callable(meta):
                meta = {}

            match_id = getattr(m, "id", None) or getattr(m, "get", lambda *_: None)("id")  # type: ignore[attr-defined]
            score = getattr(m, "score", None) or getattr(m, "get", lambda *_: 0.0)("score")  # type: ignore[attr-defined]

            vector_results.append(
                {
                    "id": match_id,
                    "score": score,
                    "metadata": meta or {},
                    "chunk_id": (meta or {}).get("chunk_id") or match_id,
                }
            )

        if bm25_results:
            fused = reciprocal_rank_fusion(
                bm25_results=bm25_results,
                vector_results=vector_results,
                alpha=hybrid_alpha,
            )

            # Keep only top_k fused results and adapt them back into a
            # Pinecone‑like structure for downstream processing.
            hybrid_matches: List[Dict[str, Any]] = []
            for r in fused[:top_k]:
                meta = r.get("metadata") or {}
                content = r.get("content")
                # Ensure we have a snippet for chunk collection / display
                if content and not meta.get("snippet"):
                    meta["snippet"] = content

                hybrid_matches.append(
                    {
                        "id": r.get("chunk_id") or r.get("id"),
                        "score": r.get("rrf_score") or r.get("bm25_score") or r.get("score", 0.0),
                        "metadata": meta,
                    }
                )

            if hybrid_matches:
                logger.info(
                    "[RAG][HYBRID] Using %d fused matches (vector + BM25)", len(hybrid_matches)
                )
                raw_matches = hybrid_matches

    # Optional cross-encoder reranking
    if use_rerank and raw_matches:
        logger.info("[RAG][RERANK] Reranking with cross-encoder")
        reranker = get_reranker()
        if rerank_model and hasattr(reranker, "model_name") and reranker.model_name != rerank_model:
            if getattr(reranker, "_model", None) is None:
                reranker.model_name = rerank_model  # type: ignore[attr-defined]
            else:
                logger.warning(
                    "[RAG][RERANK] Reranker already loaded with model %s; ignoring requested %s",
                    reranker.model_name,
                    rerank_model,
                )

        documents: List[Dict[str, Any]] = []
        for m in raw_matches:
            if isinstance(m, dict):
                meta = m.get("metadata") or {}
                chunk_id = m.get("chunk_id") or m.get("id")
                score = m.get("score", 0.0)
            else:
                meta = getattr(m, "metadata", {}) or {}
                chunk_id = getattr(m, "id", None)
                score = getattr(m, "score", 0.0)

            doc_text = (
                meta.get("chunk_text")
                or meta.get("snippet")
                or meta.get("text")
                or meta.get("content")
                or ""
            )

            documents.append(
                {
                    "id": chunk_id,
                    "chunk_text": doc_text,
                    "score": score,
                    "metadata": meta,
                }
            )

        try:
            reranked = reranker.rerank(user_query, documents, top_k=top_k)
            raw_matches = []
            for r in reranked:
                meta = r.get("metadata") or {}
                raw_matches.append(
                    {
                        "id": r.get("id"),
                        "score": r.get("rerank_score"),
                        "metadata": meta,
                    }
                )
        except Exception as e:
            logger.warning("[RAG][RERANK] Rerank failed, using original scores: %s", e)

    if not raw_matches:
        return RAGQueryResult(
            query=user_query,
            matches=[],
            filtered_matches=[],
            context_chunks=[],
            answer="No matches returned from Pinecone.",
        )

    matches = [summarize_match(m) for m in raw_matches]
    filtered = filter_matches(matches, tag=tag)

    if not filtered:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=[],
            context_chunks=[],
            answer="No matches left after filtering (tag).",
        )

    # Use hybrid chunk collection if Postgres is enabled
    if use_postgres:
        chunks = collect_hybrid_chunks(filtered, prefer_postgres=True)
    else:
        chunks = collect_chunks(filtered)
    
    if not chunks:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=filtered,
            context_chunks=[],
            answer="No Notion chunks could be retrieved for the filtered matches.",
        )

    if not generate_answer:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=filtered,
            context_chunks=chunks,
            answer="Answer generation disabled (generate_answer=False).",
        )

    metadata_list = [m.metadata for m in filtered]
    answer, citations = synthesize_answer_with_citations(user_query, chunks, metadata_list)

    # Check for hallucinations if enabled
    if check_hallucination and answer:
        logger.info("[RAG] Checking for hallucinations")
        groundedness = check_groundedness(answer, chunks)
        groundedness_score = groundedness.score
        unsupported_claims = groundedness.unsupported_claims
    else:
        groundedness_score = None
        unsupported_claims = []

    # Evaluate response if enabled
    if evaluate and answer:
        logger.info("[RAG] Running evaluation metrics")
        eval_result = evaluate_rag_response(user_query, answer, chunks)
    else:
        eval_result = None

    # Add provenance summary
    source_counts: Dict[str, int] = {}
    for m in filtered:
        src = m.metadata.get("source_type") or m.metadata.get("source") or "Unknown"
        source_counts[src] = source_counts.get(src, 0) + 1

    provenance_lines = [
        "\nSources used:",
        *[f"- {src}: {cnt} chunk(s)" for src, cnt in sorted(source_counts.items())],
    ]

    answer_with_provenance = answer + "\n\n" + "\n".join(provenance_lines)

    result = RAGQueryResult(
        query=user_query,
        matches=matches,
        filtered_matches=filtered,
        context_chunks=chunks,
        answer=answer_with_provenance,
        citations=citations,
        groundedness_score=groundedness_score,
        unsupported_claims=unsupported_claims,
        evaluation=eval_result,
    )

    # Store in cache
    if use_cache and generate_answer:
        cache = get_semantic_cache()
        cache.set(user_query, result)

    return result


def signature_similarity_query(
    dataset_page_id: str,
    top_k: int = 10,
) -> List[Dict[str, Any]]:
    """
    Given a dataset page ID (Experimental Data Assets), compute and return
    the top_k matching signatures with their scores and overlap metrics.

    Args:
        dataset_page_id: Notion page ID of dataset (with or without dashes)
        top_k: Number of top signatures to return (default: 10)

    Returns:
        List of dicts, each containing:
          - signature_page_id
          - signature_name
          - signature_short_id (if available)
          - disease (if available)
          - matrix (if available)
          - score
          - overlap_fraction
          - matched_components
          - missing_components
          - conflicting_components
    """
    logger.info(
        "[RAG][SIGNATURE-SCORE] Computing signature similarity for dataset %s (top_k=%d)",
        dataset_page_id,
        top_k,
    )

    try:
        # DEPRECATED: Notion support removed
        # TODO: Phase 3 - Update to use Postgres dataset_id instead of Notion page_id
        logger.warning(
            "[RAG][SIGNATURE-SCORE] signature_similarity_query() using Notion page_id is deprecated. "
            "Use Postgres dataset_id instead."
        )
        
        # Stub: Return empty for now
        # In Phase 3, this should:
        # 1. Accept dataset_id: UUID instead of dataset_page_id: str
        # 2. Fetch dataset from Postgres
        # 3. Extract features from Postgres or mwTab API
        # 4. Use fetch_all_signatures_from_postgres() and load_signature_from_postgres()
        return []
    except Exception as e:
        logger.error(
            "[RAG][SIGNATURE-SCORE] ERROR: Failed to compute signature similarity for dataset %s: %r",
            dataset_page_id,
            e,
        )
        return []

