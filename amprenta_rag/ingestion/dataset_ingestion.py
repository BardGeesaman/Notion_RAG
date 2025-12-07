# amprenta_rag/ingestion/dataset_ingestion.py
"""
Dataset ingestion module.

Handles ingestion of experimental datasets from Notion into Pinecone.
Supports mwTab data extraction, species extraction, signature matching,
and automatic metadata population.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_mwtab, link_features_to_notion_items)
from amprenta_rag.ingestion.metadata_semantic import \
    get_dataset_semantic_metadata
from amprenta_rag.ingestion.mwtab_extraction import (
    extract_mwtab_from_page_content, extract_metadata_from_mwtab,
    extract_study_id_from_page_properties, fetch_mwtab_from_api)
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.signature_integration import \
    detect_and_ingest_signatures_from_content
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset, map_raw_lipid_to_canonical_species,
    update_dataset_with_signature_matches)
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def ingest_dataset(page_id: str, force: bool = False) -> None:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously ingested datasets from Notion.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Remove this function or refactor to use Postgres UUIDs.
    
    Args:
        page_id: Notion page ID (ignored)
        force: If True, re-ingest (ignored)
    """
    logger.warning(
        "[INGEST][DATASET] ingest_dataset() deprecated - Notion support removed. "
        "Use Postgres-based ingestion instead."
    )
    return

    # Stub: All Notion-dependent code removed
    # Original function body removed - use Postgres-based ingestion instead
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error extracting content for %s: %r",
            page_id,
            e,
        )
        raise

    if not full_text or len(full_text.strip()) < 50:
        logger.info(
            "[INGEST][DATASET] Dataset %s has very little text; skipping.", page_id
        )
        return

    # Semantic + signature metadata
    try:
        base_meta = get_dataset_semantic_metadata(page)
        logger.info(
            "[INGEST][DATASET] Lipid signatures for %s: %r",
            page_id,
            base_meta.get("lipid_signatures"),
        )
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error reading semantic metadata for %s: %r",
            page_id,
            e,
        )
        raise

    # Chunk and embed
    chunks = chunk_text(full_text)
    if not chunks:
        logger.info("[INGEST][DATASET] No chunks produced for %s; skipping.", page_id)
        return

    logger.info(
        "[INGEST][DATASET] Generated %d chunk(s) for dataset %s",
        len(chunks),
        page_id,
    )

    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error embedding chunks for %s: %r",
            page_id,
            e,
        )
        raise

    index = get_pinecone_index()
    cfg = get_config()

    vectors: List[Dict[str, Any]] = []
    embedding_ids: List[str] = []

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id.replace('-', '')}_chunk_{order:03d}"
        embedding_ids.append(chunk_id)

        snippet = textwrap.shorten(chunk, width=300)

        meta: Dict[str, Any] = {
            **base_meta,
            "source": "Dataset",
            "source_type": "Dataset",
            "dataset_page_id": page_id.replace("-", ""),
            "title": dataset_title,
            "snippet": snippet,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(meta),
            }
        )

    if not vectors:
        logger.info("[INGEST][DATASET] No vectors to upsert for %s; skipping.", page_id)
        return

    logger.info(
        "[INGEST][DATASET] Upserting %d vectors into Pinecone for dataset %s",
        len(vectors),
        page_id,
    )

    # Batch upserts to avoid Pinecone size limits (2MB per request)
    # Rough estimate: ~15-20KB per vector (1536 dims * 4 bytes + metadata)
    # Use batches of ~100 vectors to stay well under 2MB limit
    batch_size = 100
    try:
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(vectors) + batch_size - 1) // batch_size

            logger.debug(
                "[INGEST][DATASET] Upserting batch %d/%d (%d vectors) for dataset %s",
                batch_num,
                total_batches,
                len(batch),
                page_id,
            )

            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

            if batch_num % 10 == 0 or batch_num == total_batches:
                logger.info(
                    "[INGEST][DATASET] Completed batch %d/%d for dataset %s",
                    batch_num,
                    total_batches,
                    page_id,
                )
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error upserting vectors for %s: %r",
            page_id,
            e,
        )
        raise

    # Update Notion page with embedding metadata after successful upsert
    try:
        update_dataset_embedding_metadata(
            page_id=canonical_page_id,
            embedding_ids=embedding_ids,
            embedding_count=len(embedding_ids),
        )
    except Exception as e:
        # Log but don't raise - Pinecone ingestion succeeded
        logger.warning(
            "[INGEST][DATASET] Failed to update Notion page metadata for %s (ingestion succeeded): %r",
            canonical_page_id,
            e,
        )

    # Extract mwTab data once for both metadata and feature extraction
    mwtab_data: Optional[Dict[str, Any]] = None
    try:
        mwtab_data = extract_mwtab_from_page_content(full_text)

        # Fallback: If extraction from page content failed, try fetching from MW API
        if not mwtab_data:
            logger.info(
                "[INGEST][MWTAB] Using fallback MW API fetch (page content extraction failed)"
            )

            # Try to extract STUDY_ID from multiple sources
            study_id = extract_study_id_from_page_properties(props, full_text)

            if study_id:
                logger.info(
                    "[INGEST][MWTAB] Attempting MW API fallback fetch for STUDY_ID: %s",
                    study_id,
                )
                mwtab_data = fetch_mwtab_from_api(study_id)
            else:
                logger.info(
                    "[INGEST][MWTAB] Could not extract STUDY_ID from page properties or content for MW API fallback"
                )
    except Exception as e:
        logger.warning(
            "[INGEST][DATASET] Error extracting mwTab data for %s: %r",
            page_id,
            e,
        )
        mwtab_data = None

    # Final check: if mwTab data is still None, log explicit error
    if not mwtab_data:
        logger.warning(
            "[INGEST][MWTAB] ERROR: Could not extract mwTab for dataset %s. Signature matching will be skipped.",
            page_id,
        )

    # Extract and update scientific metadata from mwTab data
    if mwtab_data:
        try:
            logger.info("[INGEST][DATASET] Extracted mwTab data for %s", page_id)
            metadata = extract_metadata_from_mwtab(mwtab_data)

            # Only update if we have at least one field to set
            has_metadata = any(
                [
                    metadata.get("model_systems"),
                    metadata.get("disease_terms"),
                    metadata.get("matrix_terms"),
                    metadata.get("methods"),
                    metadata.get("summary"),
                    metadata.get("results"),
                    metadata.get("conclusions"),
                    metadata.get("data_origin"),
                    metadata.get("dataset_source_type"),
                    metadata.get("source_url"),
                ]
            )

            if has_metadata:
                update_dataset_scientific_metadata(
                    page_id=canonical_page_id,
                    metadata=metadata,
                )
            else:
                logger.info(
                    "[INGEST][DATASET] No extractable metadata found in mwTab for %s",
                    page_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][DATASET] Error extracting/updating scientific metadata for %s: %r",
                page_id,
                e,
            )

    # Extract and link metabolite features from mwTab
    if mwtab_data:
        try:
            feature_names = extract_features_from_mwtab(mwtab_data)
            if feature_names:
                link_features_to_notion_items(
                    feature_names=feature_names,
                    item_page_id=canonical_page_id,
                    item_type="dataset",
                )
            else:
                logger.debug(
                    "[INGEST][DATASET] No metabolite features extracted from mwTab for %s",
                    page_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][DATASET] Error extracting/linking features for %s: %r",
                page_id,
                e,
            )

    # Detect and ingest signatures from dataset content
    try:
        # Combine page text and mwTab JSON as text for signature detection
        all_text_content = full_text
        if mwtab_data:
            # Convert mwTab JSON to text representation for signature detection
            mwtab_text = json.dumps(mwtab_data, indent=2)
            all_text_content = f"{full_text}\n\n{mwtab_text}"

        # Build source metadata for signature inference
        source_metadata = {
            "diseases": base_meta.get("diseases", []),
            "matrix": base_meta.get("matrix", []),
            "model_systems": base_meta.get("model_systems", []),
        }

        # No attached files for datasets typically (mwTab is in page content)
        attachment_paths: List[Path] = []

        detect_and_ingest_signatures_from_content(
            all_text_content=all_text_content,
            attachment_paths=attachment_paths,
            source_page_id=canonical_page_id,
            source_type="dataset",
            source_metadata=source_metadata,
            source_name=dataset_title,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][DATASET] Error detecting/ingesting signatures for %s: %r",
            page_id,
            e,
        )
        # Non-blocking - continue

    # Match dataset against existing signatures
    cfg = get_config()
    if cfg.pipeline.enable_signature_scoring:
        logger.info(
            "[INGEST][SIGNATURE-MATCH] Signature scoring enabled, mwTab data present: %s",
            mwtab_data is not None,
        )
    if cfg.pipeline.enable_signature_scoring and mwtab_data:
        try:
            # Extract dataset species from mwTab metabolite data
            dataset_species_set: set[str] = set()

            # Extract from mwTab MS_METABOLITE_DATA if available
            metabolite_sections = [
                "MS_METABOLITE_DATA",
                "GC_METABOLITE_DATA",
                "LC_METABOLITE_DATA",
                "METABOLITE_DATA",
            ]
            for section_key in metabolite_sections:
                if section_key in mwtab_data:
                    section = mwtab_data.get(section_key, {})
                    data_array = section.get("Data", [])
                    if isinstance(data_array, list):
                        for row in data_array:
                            if isinstance(row, dict):
                                # Look for metabolite name keys
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
                                                # Map to canonical species if possible
                                                canonical = (
                                                    map_raw_lipid_to_canonical_species(
                                                        raw_name
                                                    )
                                                )
                                                dataset_species_set.add(
                                                    canonical if canonical else raw_name
                                                )

            # Find matching signatures
            if dataset_species_set:
                logger.info(
                    "[INGEST][SIGNATURE-MATCH] Matching dataset %s against signatures (%d species)",
                    page_id,
                    len(dataset_species_set),
                )

                # Determine omics type from dataset properties
                omics_type = None
                omics_type_prop = props.get("Omics Type", {}).get("select")
                if omics_type_prop:
                    omics_type = omics_type_prop.get("name")

                matches = find_matching_signatures_for_dataset(
                    dataset_species=dataset_species_set,  # Legacy support
                    overlap_threshold=cfg.pipeline.signature_overlap_threshold,
                    dataset_page_id=canonical_page_id,  # Multi-omics support
                    omics_type=omics_type,
                )

                if matches:
                    logger.info(
                        "[INGEST][SIGNATURE-MATCH] Found %d matching signature(s) for dataset %s",
                        len(matches),
                        page_id,
                    )

                    # Update dataset page with matches
                    update_dataset_with_signature_matches(
                        dataset_page_id=canonical_page_id,
                        matches=matches,
                    )
                else:
                    logger.info(
                        "[INGEST][SIGNATURE-MATCH] No signature matches found for dataset %s",
                        page_id,
                    )
            else:
                logger.info(
                    "[INGEST][SIGNATURE-MATCH] No species extracted from dataset %s for matching",
                    page_id,
                )

        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Error matching signatures for dataset %s: %r",
                page_id,
                e,
            )
            # Non-blocking - continue

    logger.info("[INGEST][DATASET] Ingestion complete for dataset %s", page_id)
