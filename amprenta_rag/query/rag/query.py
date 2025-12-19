"""
Main RAG query orchestration.

This module provides the main query functions that orchestrate the complete
RAG pipeline: Pinecone queries, match processing, chunk collection, and answer synthesis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from amprenta_rag.query.semantic_cache import get_semantic_cache

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.rate_limit import get_rate_limiter
from amprenta_rag.query.bm25_search import bm25_search, reciprocal_rank_fusion
from amprenta_rag.query.pinecone_query import build_meta_filter, query_pinecone
from amprenta_rag.query.rag.chunk_collection import collect_chunks
from amprenta_rag.rag.hybrid_chunk_collection import collect_hybrid_chunks
from amprenta_rag.query.rag.match_processing import filter_matches, summarize_match
from amprenta_rag.query.rag.models import RAGQueryResult
from amprenta_rag.query.rag.synthesis import synthesize_answer_with_citations
from amprenta_rag.query.reranker import get_reranker
from amprenta_rag.query.hyde import generate_hypothetical_answer
from amprenta_rag.query.hallucination import check_groundedness
from amprenta_rag.query.evaluation import evaluate_rag_response
from amprenta_rag.query.agent import agentic_rag
from amprenta_rag.query.trust_scoring import weight_results_by_trust, get_trust_summary

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
    user_id: Optional[str] = None,
    model: str = "gpt-4o",
    use_trust_scoring: bool = False,
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
    # Rate limiting check (30 queries per minute)
    if user_id:
        rate_limiter = get_rate_limiter()
        if not rate_limiter.check_rate_limit(user_id, "query_rag", limit=30, window_seconds=60):
            remaining = rate_limiter.get_remaining(user_id, "query_rag", limit=30, window_seconds=60)
            error_msg = f"Rate limit exceeded. Please wait before making another query. ({remaining} requests remaining)"
            logger.warning("[RAG] Rate limit exceeded for user %s", user_id)
            return RAGQueryResult(
                query=user_query,
                answer=f"Error: {error_msg}",
                matches=[],
                filtered_matches=[],
                context_chunks=[],
                citations=[],
            )

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

    # Apply trust scoring if enabled
    if use_trust_scoring and raw_matches:
        logger.info("[RAG] Applying trust scoring to matches")
        # Convert raw_matches to dict format for trust scoring
        match_dicts = []
        for m in raw_matches:
            if isinstance(m, dict):
                match_dicts.append(m)
            else:
                # Convert object to dict
                meta = getattr(m, "metadata", None) or {}
                match_dicts.append({
                    "id": getattr(m, "id", None) or "",
                    "score": getattr(m, "score", None) or 0.0,
                    "metadata": meta,
                })

        # Apply trust weighting
        weighted_dicts = weight_results_by_trust(match_dicts)

        # Convert back to original format (preserve order)
        raw_matches = []
        for wd in weighted_dicts:
            # Find original match and update score
            original_id = wd.get("id")
            for orig_m in match_dicts:
                orig_id = orig_m.get("id") if isinstance(orig_m, dict) else getattr(orig_m, "id", None)
                if orig_id == original_id:
                    # Update score in original match
                    if isinstance(orig_m, dict):
                        orig_m["score"] = wd.get("final_score", orig_m.get("score", 0.0))
                        orig_m["trust_score"] = wd.get("trust_score", 0.0)
                    else:
                        orig_m.score = wd.get("final_score", getattr(orig_m, "score", 0.0))
                    raw_matches.append(orig_m)
                    break

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
    answer, citations = synthesize_answer_with_citations(user_query, chunks, metadata_list, model=model)

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

    # Get trust summary if trust scoring enabled
    trust_summary = None
    if use_trust_scoring and filtered:
        # Convert filtered matches to dict format for trust summary
        match_dicts = []
        for m in filtered:
            match_dicts.append({
                "score": m.score,
                "metadata": m.metadata,
            })
        trust_summary = get_trust_summary(match_dicts)

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
        trust_summary=trust_summary,
    )

    # Store in cache
    if use_cache and generate_answer:
        cache = get_semantic_cache()
        cache.set(user_query, result)

    return result


