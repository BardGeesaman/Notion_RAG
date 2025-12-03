"""
Internal polar metabolomics dataset ingestion pipeline.

This module handles ingestion of internal Amprenta metabolomics files (CSV/TSV),
including metabolite normalization, Notion integration, and RAG embedding.
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
    "ingest_metabolomics_file",
    "normalize_metabolite_name",
    "extract_metabolite_set_from_file",
]


def normalize_metabolite_name(raw: str) -> str:
    """
    Normalize metabolite names into a canonical form (best-effort).

    Phase 1: light normalization + synonym cleanup.

    Args:
        raw: Raw metabolite name from file

    Returns:
        Normalized/canonical metabolite name
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Remove adducts: [M+H], [M-H], [M+Na], etc.
    normalized = re.sub(r"\[M[+-][A-Za-z0-9]+\]", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\[M\+H\]\+", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\[M-H\]-", "", normalized, flags=re.IGNORECASE)
    # Remove standalone charge indicators: +, -, +1, -1, etc.
    normalized = re.sub(r"\s*[+-]\d*\s*$", "", normalized)
    normalized = re.sub(r"\s*\+\s*$", "", normalized)
    normalized = re.sub(r"\s*-\s*$", "", normalized)

    # Remove trailing annotations: (pos), (neg), (+), (-)
    normalized = re.sub(r"\s*\(pos\)", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\s*\(neg\)", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\s*\(\+\)", "", normalized)
    normalized = re.sub(r"\s*\(-\)", "", normalized)

    # Replace underscores with spaces
    normalized = normalized.replace("_", " ")

    # Trim again after removals
    normalized = normalized.strip()

    # Basic synonym mapping (optional, minimal for Phase 1)
    synonym_map = {
        "l-glutamic acid": "Glutamate",
        "glutamic acid": "Glutamate",
        "l-glutamine": "Glutamine",
        "glutamine": "Glutamine",
        "l-serine": "Serine",
        "serine": "Serine",
        "l-alanine": "Alanine",
        "alanine": "Alanine",
        "l-aspartic acid": "Aspartate",
        "aspartic acid": "Aspartate",
        "l-aspartate": "Aspartate",
        "aspartate": "Aspartate",
        "l-lysine": "Lysine",
        "lysine": "Lysine",
        "l-arginine": "Arginine",
        "arginine": "Arginine",
        "l-proline": "Proline",
        "proline": "Proline",
        "l-valine": "Valine",
        "valine": "Valine",
        "l-leucine": "Leucine",
        "leucine": "Leucine",
        "l-isoleucine": "Isoleucine",
        "isoleucine": "Isoleucine",
        "l-methionine": "Methionine",
        "methionine": "Methionine",
        "l-phenylalanine": "Phenylalanine",
        "phenylalanine": "Phenylalanine",
        "l-tyrosine": "Tyrosine",
        "tyrosine": "Tyrosine",
        "l-tryptophan": "Tryptophan",
        "tryptophan": "Tryptophan",
        "l-cysteine": "Cysteine",
        "cysteine": "Cysteine",
        "l-threonine": "Threonine",
        "threonine": "Threonine",
        "l-histidine": "Histidine",
        "histidine": "Histidine",
    }

    normalized_lower = normalized.lower()
    if normalized_lower in synonym_map:
        normalized = synonym_map[normalized_lower]

    # Case normalization: Title Case for common metabolites
    # Keep as-is if already looks formatted
    if normalized and not normalized[0].isupper():
        # Convert to title case
        normalized = normalized.title()

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][METABOLOMICS] Normalized metabolite '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][METABOLOMICS] Keeping metabolite name as-is: '%s'",
            original,
        )

    return normalized


def extract_metabolite_set_from_file(file_path: str) -> tuple[Set[str], int]:
    """
    Load CSV/TSV and return a set of normalized metabolite names.

    Args:
        file_path: Path to the metabolomics file

    Returns:
        Tuple of (set of normalized metabolites, total number of raw rows)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][METABOLOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][METABOLOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][METABOLOMICS] File %s is empty", file_path)
        return set(), 0

    # Detect metabolite identity column
    metabolite_column = None
    candidate_names = [
        "metabolite",
        "Metabolite",
        "compound",
        "Compound",
        "name",
        "Name",
        "molecule",
        "Molecule",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            metabolite_column = col
            break

    if not metabolite_column:
        # Try case-insensitive match
        for col in df.columns:
            if any(cand.lower() in col.lower() for cand in candidate_names):
                metabolite_column = col
                break

    if not metabolite_column:
        logger.error(
            "[INGEST][METABOLOMICS] Could not find metabolite identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find metabolite identity column in {file_path}"
        )

    logger.info(
        "[INGEST][METABOLOMICS] Using column '%s' for metabolite identity in file %s",
        metabolite_column,
        file_path,
    )

    # Extract and normalize metabolites
    metabolite_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(metabolite_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_metabolite_name(raw_name_str)
        if normalized:
            metabolite_set.add(normalized)

    logger.info(
        "[INGEST][METABOLOMICS] Extracted %d unique metabolites from %d rows in file %s",
        len(metabolite_set),
        total_rows,
        file_path,
    )

    return metabolite_set, total_rows


def create_metabolomics_dataset_page(
    file_path: str,
    metabolite_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal metabolomics dataset.

    Args:
        file_path: Path to the metabolomics file
        metabolite_count: Number of unique normalized metabolites
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Metabolomics",
        entity_count=metabolite_count,
        raw_rows=raw_rows,
        entity_name="metabolites",
    )






def embed_metabolomics_dataset(
    page_id: str,
    dataset_name: str,
    metabolites: Set[str],
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Embed metabolomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        metabolites: Set of normalized metabolites
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    cfg = get_config()

    try:
        # Build text representation
        text_parts = [
            f"Dataset: {dataset_name}",
            "Type: Metabolomics (internal)",
            "Data Origin: Internal – Amprenta",
        ]

        if program_ids:
            text_parts.append(f"Programs: {len(program_ids)} program(s) linked")

        if experiment_ids:
            text_parts.append(f"Experiments: {len(experiment_ids)} experiment(s) linked")

        text_parts.append("")
        text_parts.append(f"Metabolites ({len(metabolites)} unique):")

        # Add metabolite list (truncate if too long)
        metabolite_list = sorted(list(metabolites))
        if len(metabolite_list) > 100:
            text_parts.append(", ".join(metabolite_list[:100]))
            text_parts.append(f"... and {len(metabolite_list) - 100} more metabolites")
        else:
            text_parts.append(", ".join(metabolite_list))

        dataset_text = "\n".join(text_parts)

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][METABOLOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_metab_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "omics_type": "Metabolomics",
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
            "[INGEST][METABOLOMICS] Generated %d chunk(s)",
            len(chunks),
        )
        logger.info(
            "[INGEST][METABOLOMICS] Upserted %d vectors to Pinecone",
            len(vectors),
        )
        logger.info(
            "[INGEST][METABOLOMICS] Updated Embedding IDs and Last Embedded for dataset %s",
            page_id,
        )

    except Exception as e:
        logger.warning(
            "[INGEST][METABOLOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical


def ingest_metabolomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal polar metabolomics dataset file.

    Args:
        file_path: Path to the metabolomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][METABOLOMICS] Starting ingestion of metabolomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Metabolomics file not found: {file_path}")

    # Extract metabolites from file
    metabolite_set, raw_rows = extract_metabolite_set_from_file(file_path)

    if not metabolite_set:
        raise ValueError(f"No metabolites extracted from file: {file_path}")

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][METABOLOMICS] Using existing dataset page %s",
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
                f"Added metabolomics file {file_basename} with {len(metabolite_set)} unique metabolites."
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
                "[INGEST][METABOLOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_metabolomics_dataset_page(
            file_path=file_path,
            metabolite_count=len(metabolite_set),
            raw_rows=raw_rows,
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Metabolomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Metabolomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Link metabolites to Metabolite Features DB
    try:
        from amprenta_rag.ingestion.feature_extraction import link_feature

        logger.info(
            "[INGEST][METABOLOMICS] Linking %d metabolites to Metabolite Features DB",
            len(metabolite_set),
        )
        linked_count = 0
        for metabolite in metabolite_set:
            try:
                link_feature("metabolite", metabolite, page_id)
                linked_count += 1
            except Exception as e:
                logger.warning(
                    "[INGEST][METABOLOMICS] Error linking metabolite '%s': %r",
                    metabolite,
                    e,
                )
        logger.info(
            "[INGEST][METABOLOMICS] Linked %d/%d metabolites to Metabolite Features DB",
            linked_count,
            len(metabolite_set),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][METABOLOMICS] Feature linking skipped (error): %r",
            e,
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    embed_metabolomics_dataset(
        page_id=page_id,
        dataset_name=dataset_name,
        metabolites=metabolite_set,
        program_ids=program_ids,
        experiment_ids=experiment_ids,
    )

    logger.info(
        "[INGEST][METABOLOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

