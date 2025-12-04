"""
Main transcriptomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
transcriptomics dataset ingestion pipeline: file parsing, Notion page creation,
feature linking, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd
import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.ingestion.transcriptomics.embedding import embed_transcriptomics_dataset
from amprenta_rag.ingestion.transcriptomics.file_parsing import extract_gene_set_from_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_transcriptomics_dataset_page(
    file_path: str,
    gene_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal transcriptomics dataset.

    Args:
        file_path: Path to the transcriptomics file
        gene_count: Number of unique normalized genes
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Transcriptomics",
        entity_count=gene_count,
        raw_rows=raw_rows,
        entity_name="genes",
    )


def ingest_transcriptomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal transcriptomics dataset file.

    Args:
        file_path: Path to the transcriptomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Starting ingestion of transcriptomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Transcriptomics file not found: {file_path}")

    # Extract genes from file
    gene_set, df = extract_gene_set_from_file(file_path)

    if not gene_set:
        raise ValueError(f"No genes extracted from file: {file_path}")

    # Find gene column name for embedding
    gene_column = None
    candidate_names = ["gene", "Gene", "gene_name", "Gene Symbol", "GeneSymbol", "symbol", "Symbol"]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            gene_column = col
            break
    if not gene_column:
        # Fallback: use first column
        gene_column = df.columns[0] if len(df.columns) > 0 else "gene"

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Using existing dataset page %s",
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
                f"Added transcriptomics DGE file {file_basename} with {len(gene_set)} unique genes."
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
                "[INGEST][TRANSCRIPTOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_transcriptomics_dataset_page(
            file_path=file_path,
            gene_count=len(gene_set),
            raw_rows=len(df),
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Transcriptomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Transcriptomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Link genes to Gene Features DB
    try:
        from amprenta_rag.ingestion.feature_extraction import link_feature

        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Linking %d genes to Gene Features DB",
            len(gene_set),
        )
        linked_count = 0
        for gene in gene_set:
            try:
                link_feature("gene", gene, page_id)
                linked_count += 1
            except Exception as e:
                logger.warning(
                    "[INGEST][TRANSCRIPTOMICS] Error linking gene '%s': %r",
                    gene,
                    e,
                )
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Linked %d/%d genes to Gene Features DB",
            linked_count,
            len(gene_set),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][TRANSCRIPTOMICS] Feature linking skipped (error): %r",
            e,
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    embed_transcriptomics_dataset(
        page_id=page_id,
        dataset_name=dataset_name,
        genes=gene_set,
        df=df,
        gene_column=gene_column,
        program_ids=program_ids,
        experiment_ids=experiment_ids,
    )

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

