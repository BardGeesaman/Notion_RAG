"""
Main RAG query orchestration.

This module provides the main query functions that orchestrate the complete
RAG pipeline: Pinecone queries, match processing, chunk collection, and answer synthesis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page
from amprenta_rag.ingestion.feature_extraction import extract_features_from_mwtab
from amprenta_rag.ingestion.mwtab_extraction import (
    extract_mwtab_from_page_content,
    fetch_mwtab_from_api,
)
from amprenta_rag.ingestion.notion_pages import extract_page_content
from amprenta_rag.ingestion.signature_matching.signature_loader import (
    fetch_all_signatures_from_notion,
    load_signature_from_notion_page,
)
from amprenta_rag.ingestion.signature_matching.matching import (
    score_signature_against_dataset,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.pinecone_query import build_meta_filter, query_pinecone
from amprenta_rag.query.rag.chunk_collection import collect_chunks
from amprenta_rag.query.rag.match_processing import filter_matches, summarize_match
from amprenta_rag.query.rag.models import RAGQueryResult
from amprenta_rag.query.rag.synthesis import synthesize_answer

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
) -> RAGQueryResult:
    """
    High-level API: run a complete RAG query and get structured results.

    This function orchestrates the full RAG pipeline:
    1. Builds metadata filter from query constraints
    2. Queries Pinecone for matches (with source type filtering)
    3. Summarizes and filters matches
    4. Retrieves chunk text from Notion
    5. Synthesizes answer using OpenAI (if enabled)

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

    Returns:
        RAGQueryResult with matches, filtered matches, context chunks, and answer
    """
    logger.info("[RAG] Querying Pinecone (top_k=%d) for: %s", top_k, user_query)
    meta_filter = build_meta_filter(
        disease=disease,
        target=target,
        lipid=lipid,
        signature=signature,
    )

    raw_matches = query_pinecone(
        user_query,
        top_k=top_k,
        meta_filter=meta_filter,
        source_types=source_types,
    )
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

    answer = synthesize_answer(user_query, chunks)

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

    return RAGQueryResult(
        query=user_query,
        matches=matches,
        filtered_matches=filtered,
        context_chunks=chunks,
        answer=answer_with_provenance,
    )


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
        # Fetch dataset page
        dataset_page = fetch_dataset_page(dataset_page_id)
        if not dataset_page:
            logger.error(
                "[RAG][SIGNATURE-SCORE] Dataset page %s not found",
                dataset_page_id,
            )
            return []

        # Extract page content to get mwTab data
        page_content = extract_page_content(dataset_page_id)
        if not page_content:
            logger.warning(
                "[RAG][SIGNATURE-SCORE] Could not extract content from dataset page %s",
                dataset_page_id,
            )
            return []

        # Extract mwTab data
        mwtab_data = extract_mwtab_from_page_content(page_content)
        if not mwtab_data:
            # Try API fallback
            from amprenta_rag.ingestion.mwtab_extraction import (
                extract_study_id_from_page_properties)
            study_id = extract_study_id_from_page_properties(
                dataset_page.get("properties", {}), page_content
            )
            if study_id:
                logger.info(
                    "[RAG][SIGNATURE-SCORE] Attempting MW API fallback for STUDY_ID: %s",
                    study_id,
                )
                mwtab_data = fetch_mwtab_from_api(study_id)

        if not mwtab_data:
            logger.warning(
                "[RAG][SIGNATURE-SCORE] No mwTab data found for dataset %s",
                dataset_page_id,
            )
            return []

        # Extract species from mwTab (reuse ingestion logic)
        dataset_species_set: set[str] = set()
        feature_names = extract_features_from_mwtab(mwtab_data)

        # Also extract from MS_METABOLITE_DATA if present
        ms_data = mwtab_data.get("MS_METABOLITE_DATA", {})
        if ms_data:
            data_rows = ms_data.get("Data", [])
            if isinstance(data_rows, list):
                for row in data_rows:
                    if isinstance(row, dict):
                        for key in row.keys():
                            if key.lower() in [
                                "metabolite",
                                "metabolite_name",
                                "compound",
                                "name",
                            ]:
                                raw_name = row.get(key)
                                if raw_name and isinstance(raw_name, str):
                                    raw_name = raw_name.strip()
                                    if raw_name:
                                        from amprenta_rag.ingestion.signature_matching import (
                                            map_raw_lipid_to_canonical_species)
                                        canonical = map_raw_lipid_to_canonical_species(
                                            raw_name
                                        )
                                        dataset_species_set.add(
                                            canonical if canonical else raw_name
                                        )

        # Add feature names to species set
        for name in feature_names:
            if name:
                from amprenta_rag.ingestion.signature_matching import (
                    map_raw_lipid_to_canonical_species)
                canonical = map_raw_lipid_to_canonical_species(name)
                dataset_species_set.add(canonical if canonical else name)

        if not dataset_species_set:
            logger.warning(
                "[RAG][SIGNATURE-SCORE] No species extracted from dataset %s",
                dataset_page_id,
            )
            return []

        logger.info(
            "[RAG][SIGNATURE-SCORE] Extracted %d species from dataset %s",
            len(dataset_species_set),
            dataset_page_id,
        )

        # Fetch all signatures from Notion
        signature_pages = fetch_all_signatures_from_notion()
        if not signature_pages:
            logger.warning(
                "[RAG][SIGNATURE-SCORE] No signatures found in Notion",
            )
            return []

        # Score each signature
        results: List[Dict[str, Any]] = []
        for sig_page in signature_pages:
            try:
                signature = load_signature_from_notion_page(sig_page)
                if not signature:
                    continue

                # Score signature against dataset
                score_result = score_signature_against_dataset(
                    signature=signature,
                    dataset_species=dataset_species_set,
                    dataset_directions=None,  # Could extract from mwTab if available
                )

                # Calculate overlap fraction
                total_components = len(signature.components)
                matched_count = len(score_result.matched_species)
                overlap_fraction = (
                    matched_count / total_components if total_components > 0 else 0.0
                )

                # Get signature metadata
                props = sig_page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                signature_name = (
                    name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                )

                # Get Short ID if available
                short_id_prop = props.get("Short ID", {})
                signature_short_id = None
                if short_id_prop.get("rich_text"):
                    signature_short_id = short_id_prop["rich_text"][0].get(
                        "plain_text", ""
                    )

                # Get disease if available
                disease_prop = props.get("Disease", {})
                disease = None
                if disease_prop.get("multi_select"):
                    disease = [d.get("name", "") for d in disease_prop["multi_select"]]

                # Get matrix if available
                matrix_prop = props.get("Matrix", {})
                matrix = None
                if matrix_prop.get("multi_select"):
                    matrix = [m.get("name", "") for m in matrix_prop["multi_select"]]

                result = {
                    "signature_page_id": sig_page.get("id", ""),
                    "signature_name": signature_name,
                    "signature_short_id": signature_short_id,
                    "disease": disease,
                    "matrix": matrix,
                    "score": score_result.total_score,
                    "overlap_fraction": overlap_fraction,
                    "matched_components": list(score_result.matched_species),
                    "missing_components": list(score_result.missing_species),
                    "conflicting_components": list(score_result.conflicting_species),
                }
                results.append(result)

                logger.debug(
                    "[RAG][SIGNATURE-SCORE] Match: %s (score=%.3f, overlap=%.2f, matched=%d, missing=%d)",
                    signature_name,
                    score_result.total_score,
                    overlap_fraction,
                    len(score_result.matched_species),
                    len(score_result.missing_species),
                )

            except Exception as e:
                logger.warning(
                    "[RAG][SIGNATURE-SCORE] Error processing signature %s: %r",
                    sig_page.get("id", ""),
                    e,
                )
                continue

        # Sort by score (descending), then by overlap_fraction (descending)
        results.sort(key=lambda x: (x["score"], x["overlap_fraction"]), reverse=True)

        # Truncate to top_k
        top_results = results[:top_k]

        logger.info(
            "[RAG][SIGNATURE-SCORE] Found %d matching signatures (returning top %d)",
            len(results),
            len(top_results),
        )

        return top_results

    except Exception as e:
        logger.error(
            "[RAG][SIGNATURE-SCORE] ERROR: Failed to compute signature similarity for dataset %s: %r",
            dataset_page_id,
            e,
        )
        return []

