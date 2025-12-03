"""
Internal proteomics dataset ingestion pipeline.

This module handles ingestion of internal Amprenta proteomics files (CSV/TSV),
including protein/gene identifier normalization, Notion integration, and RAG embedding.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import (
    fetch_dataset_page, update_dataset_embedding_metadata)
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page, create_omics_dataset_page,
    link_to_programs_and_experiments)
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "ingest_proteomics_file",
    "normalize_protein_identifier",
    "extract_protein_set_from_file",
]


def normalize_protein_identifier(raw: str) -> str:
    """
    Normalize protein or gene identifiers into a canonical symbol-like format.

    Phase 1: basic mapping + cleaning.

    Handles:
    - FASTA prefixes: sp|P12345|TP53_HUMAN → TP53
    - Isoform suffixes: Q9Y6K9-2 → Q9Y6K9
    - Gene symbols: actb → ACTB
    - Annotations: TP53 (Human) → TP53

    Args:
        raw: Raw protein/gene identifier from file

    Returns:
        Normalized/canonical identifier
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Remove FASTA prefixes: sp|P12345|TP53_HUMAN → TP53_HUMAN
    # Pattern: sp|UNIPROT_ID|GENE_SPECIES or tr|...|...
    fasta_match = re.match(r"^(sp|tr|ref)\|[^|]+\|([^|]+)", normalized)
    if fasta_match:
        normalized = fasta_match.group(2)
        logger.debug(
            "[INGEST][PROTEOMICS] Extracted from FASTA format: '%s' -> '%s'",
            original,
            normalized,
        )

    # Extract gene symbol from format like TP53_HUMAN → TP53
    # Or P04637_TP53_HUMAN → TP53
    if "_" in normalized:
        parts = normalized.split("_")
        # If last part looks like species (HUMAN, MOUSE, etc.), remove it
        if len(parts) > 1 and parts[-1].upper() in [
            "HUMAN",
            "MOUSE",
            "RAT",
            "BOVIN",
            "PIG",
            "CHICK",
        ]:
            normalized = "_".join(parts[:-1])
        # If we have something like P04637_TP53, try to extract TP53
        if len(parts) > 1:
            # Check if first part looks like UniProt ID (alphanumeric, 6-10 chars)
            if re.match(r"^[A-Z0-9]{6,10}$", parts[0].upper()):
                # Use the second part as gene symbol
                if len(parts) > 1:
                    normalized = parts[1]

    # Remove isoform suffixes: Q9Y6K9-2 → Q9Y6K9
    # But keep if it's part of a gene symbol pattern
    if re.match(r"^[A-Z0-9]+-\d+$", normalized.upper()):
        # Looks like UniProt ID with isoform, remove suffix
        normalized = re.sub(r"-\d+$", "", normalized)
    elif "-" in normalized and not re.match(r"^[A-Z]+\d+[A-Z]*$", normalized.upper()):
        # Has hyphen but not a simple gene symbol, might be isoform
        # Only remove if it's at the end and looks like a number
        normalized = re.sub(r"-\d+$", "", normalized)

    # Remove bracketed annotations: TP53 (Human) → TP53
    normalized = re.sub(r"\s*\([^)]+\)", "", normalized)
    normalized = re.sub(r"\s*\[[^\]]+\]", "", normalized)

    # Convert to uppercase for gene symbols (if it looks like a gene symbol)
    # Gene symbols are typically 1-10 uppercase letters/numbers
    if re.match(r"^[A-Za-z0-9]{1,10}$", normalized):
        normalized = normalized.upper()
    # If it's a UniProt ID (6-10 alphanumeric), keep uppercase
    elif re.match(r"^[A-Z0-9]{6,10}$", normalized.upper()):
        normalized = normalized.upper()

    # Replace underscores with hyphens only if it looks gene-like
    # (Not for UniProt IDs which may have underscores)
    if re.match(r"^[A-Z]+\d*[A-Z]*$", normalized.upper()):
        normalized = normalized.replace("_", "-")

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][PROTEOMICS] Normalized '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][PROTEOMICS] Keeping identifier as-is: '%s'",
            original,
        )

    return normalized


def extract_protein_set_from_file(file_path: str) -> tuple[Set[str], int]:
    """
    Load CSV/TSV and return a set of normalized protein/gene identifiers.

    Args:
        file_path: Path to the proteomics file

    Returns:
        Tuple of (set of normalized proteins, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][PROTEOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][PROTEOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][PROTEOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect protein/gene identity column
    protein_column = None
    candidate_names = [
        "Protein",
        "protein",
        "Protein ID",
        "ProteinID",
        "protein_id",
        "Gene",
        "gene",
        "Gene Symbol",
        "GeneSymbol",
        "gene_name",
        "protein_name",
        "Protein Name",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            protein_column = col
            break

    if not protein_column:
        # Try case-insensitive match
        for col in df.columns:
            col_lower = col.lower()
            if any(
                cand.lower() in col_lower or col_lower in cand.lower()
                for cand in candidate_names
            ):
                protein_column = col
                break

    if not protein_column:
        logger.error(
            "[INGEST][PROTEOMICS] Could not find protein/gene identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find protein/gene identity column in {file_path}"
        )

    logger.info(
        "[INGEST][PROTEOMICS] Using column '%s' for protein identity in file %s",
        protein_column,
        file_path,
    )

    # Extract and normalize proteins
    protein_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(protein_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_protein_identifier(raw_name_str)
        if normalized:
            protein_set.add(normalized)

    logger.info(
        "[INGEST][PROTEOMICS] Extracted %d unique proteins from %d rows in file %s",
        len(protein_set),
        total_rows,
        file_path,
    )

    return protein_set, total_rows


def create_proteomics_dataset_page(
    file_path: str,
    protein_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal proteomics dataset.

    Args:
        file_path: Path to the proteomics file
        protein_count: Number of unique normalized proteins
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Proteomics",
        entity_count=protein_count,
        raw_rows=raw_rows,
        entity_name="proteins",
    )




def embed_proteomics_dataset(
    page_id: str,
    dataset_name: str,
    proteins: Set[str],
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Embed proteomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        proteins: Set of normalized proteins
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    cfg = get_config()

    try:
        # Build text representation
        text_parts = [
            f"Dataset: {dataset_name}",
            "Type: Proteomics (internal)",
            "Data Origin: Internal – Amprenta",
        ]

        if program_ids:
            text_parts.append(f"Programs: {len(program_ids)} program(s) linked")

        if experiment_ids:
            text_parts.append(f"Experiments: {len(experiment_ids)} experiment(s) linked")

        text_parts.append("")
        text_parts.append(f"Proteins ({len(proteins)} unique):")

        # Add protein list (truncate if too long)
        protein_list = sorted(list(proteins))
        if len(protein_list) > 100:
            text_parts.append(", ".join(protein_list[:100]))
            text_parts.append(f"... and {len(protein_list) - 100} more proteins")
        else:
            text_parts.append(", ".join(protein_list))

        dataset_text = "\n".join(text_parts)

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][PROTEOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_prot_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "omics_type": "Proteomics",
                "dataset_page_id": page_id.replace("-", ""),
                "dataset_name": dataset_name,
                "data_origin": "Internal – Amprenta",
                "snippet": chunk[:300],
            }

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        # Batch upsert
        batch_size = 100
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

        # Update Notion with embedding metadata
        embedding_ids = [v["id"] for v in vectors]
        update_dataset_embedding_metadata(
            page_id=page_id,
            embedding_ids=embedding_ids,
            embedding_count=len(embedding_ids),
        )

        logger.info(
            "[INGEST][PROTEOMICS] Generated %d chunk(s)",
            len(chunks),
        )
        logger.info(
            "[INGEST][PROTEOMICS] Upserted %d vectors to Pinecone",
            len(vectors),
        )
        logger.info(
            "[INGEST][PROTEOMICS] Updated Embedding IDs and Last Embedded for dataset %s",
            page_id,
        )

    except Exception as e:
        logger.warning(
            "[INGEST][PROTEOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical


def ingest_proteomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal proteomics dataset file.

    Args:
        file_path: Path to the proteomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][PROTEOMICS] Starting ingestion of proteomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Proteomics file not found: {file_path}")

    # Extract proteins from file
    protein_set, raw_rows = extract_protein_set_from_file(file_path)

    if not protein_set:
        raise ValueError(f"No proteins extracted from file: {file_path}")

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][PROTEOMICS] Using existing dataset page %s",
            page_id,
        )
        # Update summary
        try:
            page = fetch_dataset_page(page_id)
            page_props = page.get("properties", {}) or {}
            summary_prop = page_props.get("Summary", {}) or {}
            existing_summary = "".join(
                rt.get("plain_text", "")
                for rt in summary_prop.get("rich_text", [])
            )
            file_basename = Path(file_path).name
            new_note = (
                f"\n\n---\n\n"
                f"Added proteomics file {file_basename} with {len(protein_set)} unique proteins."
            )
            combined_summary = existing_summary + new_note

            # Update page
            url = f"{get_config().notion.base_url}/pages/{page_id}"
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json={"properties": {"Summary": {"rich_text": [{"text": {"content": combined_summary}}]}}},
                timeout=30,
            )
            resp.raise_for_status()
        except Exception as e:
            logger.warning(
                "[INGEST][PROTEOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_proteomics_dataset_page(
            file_path=file_path,
            protein_count=len(protein_set),
            raw_rows=raw_rows,
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Proteomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Proteomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Optional: Protein Features DB integration (Phase 2)
    try:
        from amprenta_rag.ingestion.feature_extraction import (
            link_features_to_notion_items)
        cfg = get_config()
        if hasattr(cfg.notion, "protein_features_db_id") and cfg.notion.protein_features_db_id:
            # Link proteins to Protein Features DB
            protein_list = list(protein_set)
            link_features_to_notion_items(
                feature_names=protein_list,
                item_page_id=page_id,
                item_type="dataset",
            )
            logger.info(
                "[INGEST][PROTEOMICS] Linked %d proteins to Protein Features DB",
                len(protein_list),
            )
    except Exception as e:
        logger.warning(
            "[INGEST][PROTEOMICS] Protein Features DB integration skipped (not configured or error): %r",
            e,
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    embed_proteomics_dataset(
        page_id=page_id,
        dataset_name=dataset_name,
        proteins=protein_set,
        program_ids=program_ids,
        experiment_ids=experiment_ids,
    )

    logger.info(
        "[INGEST][PROTEOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

