"""
Postgres-based signature detection and ingestion.

Detects and creates signatures from content using Postgres instead of Notion.
"""

from __future__ import annotations

import tempfile
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.ingestion.signature_detection import (
    detect_signature_keywords,
    extract_embedded_signature_table,
    extract_signature_from_text_table,
    find_attached_signature_files,
    infer_signature_metadata_from_source,
    save_extracted_signature_to_file,
)
from amprenta_rag.ingestion.postgres_signature_creation import (
    create_signature_from_file_in_postgres,
    link_signature_to_postgres_source,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import load_signature_from_tsv

logger = get_logger(__name__)


def detect_and_ingest_signatures_from_postgres_content(
    all_text_content: str,
    attachment_paths: List[Path],
    source_type: str,
    source_metadata: Dict[str, Any],
    source_name: Optional[str] = None,
    source_dataset_id: Optional[UUID] = None,
    source_experiment_id: Optional[UUID] = None,
) -> Dict[str, Any]:
    """
    Detect and ingest signatures from source content using Postgres.
    
    This is the Postgres-first version that works without Notion page IDs.
    
    Args:
        all_text_content: Combined text content from the source
        attachment_paths: List of paths to attached files
        source_type: One of "literature", "dataset", "email", "experiment"
        source_metadata: Metadata dict from the source
        source_name: Optional name for the source (for signature naming)
        source_dataset_id: Optional Postgres UUID of dataset (for linking)
        source_experiment_id: Optional Postgres UUID of experiment (for linking)
        
    Returns:
        Dictionary with counts: {"detected": N, "ingested": M, "errors": K, "signature_ids": [...]}
    """
    counts: Dict[str, Any] = {"detected": 0, "ingested": 0, "errors": 0, "signature_ids": []}
    
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
                "[POSTGRES-SIGNATURE-DETECT] Signature keywords detected in %s source",
                source_type,
            )
            counts["detected"] += 1
            
            # Try to extract embedded table
            embedded_table = extract_embedded_signature_table(all_text_content)
            if embedded_table:
                signature_candidates.append(("embedded_table", embedded_table))
                logger.debug(
                    "[POSTGRES-SIGNATURE-DETECT] Found embedded signature table",
                )
            
            # Try to extract from text table
            text_table = extract_signature_from_text_table(all_text_content)
            if text_table and text_table.get("components"):
                signature_candidates.append(("text_table", text_table))
                logger.debug(
                    "[POSTGRES-SIGNATURE-DETECT] Found text table signature",
                )
        
        # 2. Check for attached signature files
        if attachment_paths:
            attached_signature_files = find_attached_signature_files(attachment_paths)
            for sig_file in attached_signature_files:
                signature_candidates.append(("file", sig_file))
                counts["detected"] += 1
                logger.debug(
                    "[POSTGRES-SIGNATURE-DETECT] Found signature file: %s",
                    sig_file.name,
                )
        
        # 3. Process each signature candidate
        for candidate_type, candidate_data in signature_candidates:
            try:
                sig_file_path: Optional[Path] = None
                
                if candidate_type == "embedded_table":
                    # Save embedded table to temporary file
                    components = candidate_data  # List[Dict[str, str]]
                    if not components or len(components) < 2:
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        sig_name = source_name or f"{source_type}_signature"
                        sig_file_path = save_extracted_signature_to_file(
                            components,
                            Path(tmpdir),
                            f"{sig_name}_signature",
                        )
                        
                elif candidate_type == "text_table":
                    # Process text table
                    components = candidate_data.get("components", [])
                    if not components or len(components) < 2:
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        sig_name = source_name or f"{source_type}_signature"
                        sig_file_path = save_extracted_signature_to_file(
                            components,
                            Path(tmpdir),
                            f"{sig_name}_signature",
                        )
                        
                elif candidate_type == "file":
                    # Direct file ingestion
                    sig_file_path = candidate_data  # Path
                
                if not sig_file_path or not sig_file_path.exists():
                    continue
                
                # Load signature from file
                signature = load_signature_from_tsv(sig_file_path)
                
                # Create signature in Postgres
                signature_model = create_signature_from_file_in_postgres(
                    signature=signature,
                    signature_type=inferred_metadata.get("signature_type", "Literature-derived"),
                    data_ownership="Public",
                    description=f"Signature extracted from {source_type}",
                    biomarker_roles=inferred_metadata.get("biomarker_roles"),
                    phenotype_axes=inferred_metadata.get("phenotype_axes"),
                )
                
                counts["ingested"] += 1
                counts["signature_ids"].append(str(signature_model.id))
                
                # Link signature to source
                if source_dataset_id:
                    link_signature_to_postgres_source(
                        signature_id=signature_model.id,
                        source_type="dataset",
                        source_id=source_dataset_id,
                    )
                elif source_experiment_id:
                    link_signature_to_postgres_source(
                        signature_id=signature_model.id,
                        source_type="experiment",
                        source_id=source_experiment_id,
                    )
                
                # Embed signature into Pinecone
                try:
                    embed_signature_with_postgres_id(
                        signature_id=signature_model.id,
                        signature=signature,
                    )
                except Exception as e:
                    logger.warning(
                        "[POSTGRES-SIGNATURE-DETECT] Error embedding signature %s: %r",
                        signature_model.id,
                        e,
                    )
                
                logger.info(
                    "[POSTGRES-SIGNATURE-DETECT] Created signature %s (ID: %s) from %s",
                    signature_model.name,
                    signature_model.id,
                    source_type,
                )
                
            except Exception as e:
                counts["errors"] += 1
                logger.warning(
                    "[POSTGRES-SIGNATURE-DETECT] Error processing signature candidate: %r",
                    e,
                )
                continue
        
    except Exception as e:
        logger.warning(
            "[POSTGRES-SIGNATURE-DETECT] Error in signature detection: %r",
            e,
        )
        counts["errors"] += 1
    
    if counts["ingested"] > 0:
        logger.info(
            "[POSTGRES-SIGNATURE-DETECT] Processed signatures: %d detected, %d ingested, %d errors",
            counts["detected"],
            counts["ingested"],
            counts["errors"],
        )
    
    return counts


def embed_signature_with_postgres_id(
    signature_id: UUID,
    signature: Signature,
) -> None:
    """
    Embed a signature into Pinecone using Postgres signature_id.
    
    Args:
        signature_id: Postgres UUID of the signature
        signature: Signature object with components
    """
    from amprenta_rag.clients.pinecone_client import get_pinecone_index
    from amprenta_rag.config import get_config
    from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
    from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
    
    cfg = get_config()
    
    try:
        # Build text representation
        signature_type = "Multi-Omics Signature" if (signature.modalities and len(signature.modalities) > 1) else "Signature"
        text_parts = [
            f"{signature_type}: {signature.name}",
        ]
        
        if signature.modalities:
            modalities_str = ", ".join(mod.title() for mod in signature.modalities)
            text_parts.append(f"Modalities: {modalities_str}")
        
        if signature.description:
            text_parts.append(f"\nDescription: {signature.description}")
        
        text_parts.append("\nComponents:")
        for comp in signature.components:
            feature_type = getattr(comp, "feature_type", "lipid")
            feature_name = getattr(comp, "feature_name", comp.species)
            direction = getattr(comp, "direction", "")
            weight = getattr(comp, "weight", 1.0)
            
            comp_text = f"- {feature_name} ({feature_type})"
            if direction:
                comp_text += f" {direction}"
            if weight and weight != 1.0:
                comp_text += f" (weight: {weight})"
            
            text_parts.append(comp_text)
        
        full_text = "\n".join(text_parts)
        
        # Chunk and embed
        chunks = chunk_text(full_text)
        if not chunks:
            logger.warning(
                "[POSTGRES-SIGNATURE-EMBED] No chunks generated for signature %s",
                signature_id,
            )
            return
        
        embeddings = embed_texts(chunks)
        
        # Prepare vectors for Pinecone
        index = get_pinecone_index()
        vectors = []
        signature_id_str = str(signature_id).replace("-", "")
        
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{signature_id_str}_chunk_{order:03d}"
            
            metadata: Dict[str, Any] = {
                "source": "Signature",
                "source_type": "Signature",
                "signature_id": str(signature_id),
                "title": signature.name,
                "snippet": textwrap.shorten(chunk, width=300),
            }
            
            if signature.modalities:
                metadata["modalities"] = signature.modalities
            
            vectors.append({
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(metadata),
            })
        
        # Upsert to Pinecone
        if vectors:
            index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
            logger.info(
                "[POSTGRES-SIGNATURE-EMBED] Embedded signature %s with %d chunk(s)",
                signature_id,
                len(vectors),
            )
        
    except Exception as e:
        logger.error(
            "[POSTGRES-SIGNATURE-EMBED] Error embedding signature %s: %r",
            signature_id,
            e,
        )
        raise

