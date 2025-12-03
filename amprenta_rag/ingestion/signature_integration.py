# amprenta_rag/ingestion/signature_integration.py

"""
Signature detection and ingestion integration for all ingestion pipelines.

This module provides a unified interface for detecting and ingesting signatures
from any source type, making it easy to integrate signature detection into
all ingestion pipelines.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional

from amprenta_rag.ingestion.signature_detection import (
    detect_signature_keywords,
    find_attached_signature_files,
    extract_embedded_signature_table,
    extract_signature_from_text_table,
    save_extracted_signature_to_file,
    infer_signature_metadata_from_source,
)
from amprenta_rag.ingestion.signature_ingestion import (
    ingest_signature_from_file,
    link_signature_to_source,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def detect_and_ingest_signatures_from_content(
    all_text_content: str,
    attachment_paths: List[Path],
    source_page_id: str,
    source_type: str,
    source_metadata: Dict[str, Any],
    source_name: Optional[str] = None,
) -> Dict[str, int]:
    """
    Detect and ingest signatures from source content.
    
    This is the main integration function that should be called from
    ingestion pipelines after Pinecone upsert but before Notion updates.
    
    Args:
        all_text_content: Combined text content from the source
        attachment_paths: List of paths to attached files
        source_page_id: Notion page ID of the source (with dashes)
        source_type: One of "literature", "dataset", "email", "experiment"
        source_metadata: Metadata dict from the source page
        source_name: Optional name for the source (for signature naming)
        
    Returns:
        Dictionary with counts: {"detected": N, "ingested": M, "errors": K}
    """
    counts = {"detected": 0, "ingested": 0, "errors": 0}
    
    if not all_text_content and not attachment_paths:
        return counts
    
    # Infer signature metadata from source
    inferred_metadata = infer_signature_metadata_from_source(
        source_type=source_type,
        source_metadata=source_metadata,
    )
    
    signature_candidates: List[tuple] = []
    
    try:
        # 1. Check for signature keywords in text
        has_keywords = detect_signature_keywords(all_text_content)
        
        if has_keywords:
            logger.debug(
                "[INGEST][SIGNATURES] Signature keywords detected in %s source %s",
                source_type,
                source_page_id,
            )
            counts["detected"] += 1
            
            # Try to extract embedded table
            embedded_table = extract_embedded_signature_table(all_text_content)
            if embedded_table:
                signature_candidates.append(("embedded_table", embedded_table))
                logger.debug(
                    "[INGEST][SIGNATURES] Found embedded signature table in %s",
                    source_page_id,
                )
            
            # Try to extract from text table
            text_table = extract_signature_from_text_table(all_text_content)
            if text_table and text_table.get("components"):
                signature_candidates.append(("text_table", text_table))
                logger.debug(
                    "[INGEST][SIGNATURES] Found text table signature in %s",
                    source_page_id,
                )
        
        # 2. Check for attached signature files
        if attachment_paths:
            attached_signature_files = find_attached_signature_files(attachment_paths)
            for sig_file in attached_signature_files:
                signature_candidates.append(("file", sig_file))
                counts["detected"] += 1
                logger.debug(
                    "[INGEST][SIGNATURES] Found signature file: %s",
                    sig_file.name,
                )
        
        # 3. Process each signature candidate
        for candidate_type, candidate_data in signature_candidates:
            try:
                if candidate_type == "embedded_table":
                    # Save embedded table to temporary file
                    components = candidate_data  # List[Dict[str, str]]
                    if not components or len(components) < 2:
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        sig_name = source_name or f"{source_type}_{source_page_id[:8]}"
                        sig_file = save_extracted_signature_to_file(
                            components,
                            Path(tmpdir),
                            f"{sig_name}_signature",
                        )
                        
                        if sig_file and sig_file.exists():
                            # Ingest signature
                            result = ingest_signature_from_file(
                                sig_file,
                                signature_type=inferred_metadata.get("signature_type", "Literature-derived"),
                                description=f"Signature extracted from {source_type}",
                                disease_context=inferred_metadata.get("disease_context"),
                                matrix=inferred_metadata.get("matrix"),
                            )
                            
                            if result.get("signature_page_id"):
                                # Link back to source
                                link_signature_to_source(
                                    result["signature_page_id"],
                                    source_page_id,
                                    source_type,
                                )
                                counts["ingested"] += 1
                                logger.info(
                                    "[INGEST][SIGNATURES] Ingested embedded signature from %s source %s",
                                    source_type,
                                    source_page_id,
                                )
                
                elif candidate_type == "text_table":
                    # Process text table
                    components = candidate_data.get("components", [])
                    if not components or len(components) < 2:
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        sig_name = source_name or f"{source_type}_{source_page_id[:8]}"
                        sig_file = save_extracted_signature_to_file(
                            components,
                            Path(tmpdir),
                            f"{sig_name}_signature",
                        )
                        
                        if sig_file and sig_file.exists():
                            result = ingest_signature_from_file(
                                sig_file,
                                signature_type=inferred_metadata.get("signature_type", "Literature-derived"),
                                description=f"Signature extracted from {source_type}",
                                disease_context=inferred_metadata.get("disease_context"),
                                matrix=inferred_metadata.get("matrix"),
                            )
                            
                            if result.get("signature_page_id"):
                                link_signature_to_source(
                                    result["signature_page_id"],
                                    source_page_id,
                                    source_type,
                                )
                                counts["ingested"] += 1
                                logger.info(
                                    "[INGEST][SIGNATURES] Ingested text table signature from %s source %s",
                                    source_type,
                                    source_page_id,
                                )
                
                elif candidate_type == "file":
                    # Direct file ingestion
                    sig_file = candidate_data  # Path
                    if not sig_file.exists():
                        continue
                    
                    result = ingest_signature_from_file(
                        sig_file,
                        signature_type=inferred_metadata.get("signature_type", "Literature-derived"),
                        description=f"Signature from {source_type} attachment",
                        disease_context=inferred_metadata.get("disease_context"),
                        matrix=inferred_metadata.get("matrix"),
                    )
                    
                    if result.get("signature_page_id"):
                        link_signature_to_source(
                            result["signature_page_id"],
                            source_page_id,
                            source_type,
                        )
                        counts["ingested"] += 1
                        logger.info(
                            "[INGEST][SIGNATURES] Ingested signature file from %s source %s: %s",
                            source_type,
                            source_page_id,
                            sig_file.name,
                        )
            
            except Exception as e:
                counts["errors"] += 1
                logger.warning(
                    "[INGEST][SIGNATURES] Error processing signature candidate in %s source %s: %r",
                    source_type,
                    source_page_id,
                    e,
                )
                # Continue with next candidate
    
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error in signature detection for %s source %s: %r",
            source_type,
            source_page_id,
            e,
        )
        counts["errors"] += 1
    
    if counts["ingested"] > 0:
        logger.info(
            "[INGEST][SIGNATURES] Processed signatures for %s source %s: "
            "%d detected, %d ingested, %d errors",
            source_type,
            source_page_id,
            counts["detected"],
            counts["ingested"],
            counts["errors"],
        )
    
    return counts

