"""
Internal transcriptomics (RNA expression) dataset ingestion pipeline.

This module handles ingestion of internal Amprenta transcriptomics files (CSV/TSV),
including gene identifier normalization, DGE data processing, Notion integration, and RAG embedding.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

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
    "ingest_transcriptomics_file",
    "normalize_gene_identifier",
    "extract_gene_set_from_file",
]


def normalize_gene_identifier(raw: str) -> str:
    """
    Normalize gene identifiers into a canonical gene symbol-like format.

    Phase 1: best-effort cleaning to HGNC-style symbols where possible.

    Handles:
    - Species suffixes: TP53_HUMAN → TP53
    - Case normalization: tp53 → TP53
    - Annotations: TP53 (Human) → TP53
    - Ensembl IDs: Keep as-is but clean whitespace

    Args:
        raw: Raw gene identifier from file

    Returns:
        Normalized/canonical gene identifier
    """
    if not raw or not isinstance(raw, str):
        return raw if raw else ""

    # Trim whitespace
    normalized = raw.strip()
    if not normalized:
        return ""

    original = normalized

    # Check if it's an Ensembl ID (ENSG0000... or similar)
    ensembl_pattern = r"^ENS[A-Z]*G\d+"
    is_ensembl = bool(re.match(ensembl_pattern, normalized, re.IGNORECASE))
    
    if is_ensembl:
        # Keep Ensembl IDs as-is but ensure uppercase and clean whitespace
        normalized = normalized.upper().strip()
        if normalized != original:
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Normalized Ensembl ID '%s' -> '%s'",
                original,
                normalized,
            )
        else:
            logger.debug(
                "[INGEST][TRANSCRIPTOMICS] Keeping Ensembl ID as-is: '%s'",
                original,
            )
        return normalized

    # Remove species suffixes: TP53_HUMAN, Tp53_mouse → TP53
    # Common species suffixes: HUMAN, MOUSE, RAT, BOVIN, PIG, CHICK
    if "_" in normalized:
        parts = normalized.split("_")
        if len(parts) > 1:
            last_part = parts[-1].upper()
            if last_part in ["HUMAN", "MOUSE", "RAT", "BOVIN", "PIG", "CHICK", "HSA", "MMU"]:
                normalized = "_".join(parts[:-1])

    # Remove bracketed annotations: TP53 (Human) → TP53
    normalized = re.sub(r"\s*\([^)]+\)", "", normalized)
    normalized = re.sub(r"\s*\[[^\]]+\]", "", normalized)

    # Remove trailing version-like suffixes if clearly non-gene: TP53.1 → TP53
    # But be careful not to remove legitimate gene symbols that end in numbers
    # Only remove if it looks like a version suffix (e.g., .1, .2, .v1)
    if re.match(r"^[A-Z]+\d*\.\d+$", normalized, re.IGNORECASE):
        normalized = re.sub(r"\.\d+$", "", normalized)

    # Convert to uppercase for gene symbols (tp53 → TP53)
    # Gene symbols are typically uppercase letters and numbers
    if re.match(r"^[A-Za-z]+\d*[A-Za-z]*$", normalized):
        normalized = normalized.upper()

    # Final trim
    normalized = normalized.strip()

    if normalized != original:
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Normalized '%s' -> '%s'",
            original,
            normalized,
        )
    else:
        logger.debug(
            "[INGEST][TRANSCRIPTOMICS] Keeping gene identifier as-is: '%s'",
            original,
        )

    return normalized


def extract_gene_set_from_file(file_path: str) -> Tuple[Set[str], pd.DataFrame]:
    """
    Load the DGE table and return a set of normalized gene identifiers and the full DataFrame.

    Args:
        file_path: Path to the transcriptomics file

    Returns:
        Tuple of (set of normalized genes, full DataFrame)
    """
    try:
        # Try to detect delimiter
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
    except Exception as e:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Error reading file %s: %r",
            file_path,
            e,
        )
        raise

    try:
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8")
    except Exception as e:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Error parsing CSV/TSV file %s: %r",
            file_path,
            e,
        )
        raise

    if df.empty:
        logger.warning("[INGEST][TRANSCRIPTOMICS] File %s is empty", file_path)
        return set(), df

    # Detect gene identity column
    gene_column = None
    candidate_names = [
        "gene",
        "Gene",
        "gene_name",
        "Gene Symbol",
        "GeneSymbol",
        "symbol",
        "Symbol",
    ]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            gene_column = col
            break

    if not gene_column:
        # Try case-insensitive match
        for col in df.columns:
            col_lower = col.lower()
            if any(
                cand.lower() in col_lower or col_lower in cand.lower()
                for cand in candidate_names
            ):
                gene_column = col
                break

    if not gene_column:
        logger.error(
            "[INGEST][TRANSCRIPTOMICS] Could not find gene identity column in file %s. "
            "Available columns: %s",
            file_path,
            list(df.columns),
        )
        raise ValueError(
            f"Could not find gene identity column in {file_path}"
        )

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Using column '%s' for gene identity in file %s",
        gene_column,
        file_path,
    )

    # Extract and normalize genes
    gene_set: Set[str] = set()
    total_rows = len(df)

    for idx, row in df.iterrows():
        raw_name = row.get(gene_column)
        if pd.isna(raw_name):
            continue

        raw_name_str = str(raw_name).strip()
        if not raw_name_str:
            continue

        normalized = normalize_gene_identifier(raw_name_str)
        if normalized:
            gene_set.add(normalized)

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Extracted %d unique genes from %d rows in file %s",
        len(gene_set),
        total_rows,
        file_path,
    )

    return gene_set, df


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




def build_dge_text_representation(
    dataset_name: str,
    genes: Set[str],
    df: pd.DataFrame,
    gene_column: str,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Build a text representation of the DGE dataset for RAG.

    Includes summary statistics and top genes by |log2FC| or p-value.

    Args:
        dataset_name: Dataset name
        genes: Set of normalized genes
        df: Full DataFrame with DGE data
        gene_column: Name of the gene column
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs

    Returns:
        Text representation of the dataset
    """
    text_parts = [
        f"Dataset: {dataset_name}",
        "Type: Transcriptomics (internal)",
        "Data Origin: Internal – Amprenta",
    ]

    if program_ids:
        text_parts.append(f"Programs: {len(program_ids)} program(s) linked")

    if experiment_ids:
        text_parts.append(f"Experiments: {len(experiment_ids)} experiment(s) linked")

    text_parts.append("")
    text_parts.append("Differential Expression Summary:")
    text_parts.append(f"- Total genes: {len(genes)}")

    # Try to find log2FC and p-value columns
    log2fc_col = None
    pval_col = None

    for col in df.columns:
        col_lower = col.lower()
        if col_lower in ["log2fc", "logfc", "log2_fc", "lfc"]:
            log2fc_col = col
        elif col_lower in ["pvalue", "p_value", "pval", "padj", "fdr", "adjpval"]:
            pval_col = col

    # Select top genes for summary (up to 50)
    if log2fc_col and pval_col:
        # Sort by absolute log2FC, then by p-value
        df_sorted = df.copy()
        df_sorted["abs_log2fc"] = df_sorted[log2fc_col].abs()
        df_sorted = df_sorted.sort_values(
            by=["abs_log2fc", pval_col], ascending=[False, True]
        )
        top_genes = df_sorted.head(50)
    elif log2fc_col:
        # Sort by absolute log2FC
        df_sorted = df.copy()
        df_sorted["abs_log2fc"] = df_sorted[log2fc_col].abs()
        df_sorted = df_sorted.sort_values(by="abs_log2fc", ascending=False)
        top_genes = df_sorted.head(50)
    elif pval_col:
        # Sort by p-value
        df_sorted = df.sort_values(by=pval_col, ascending=True)
        top_genes = df_sorted.head(50)
    else:
        # No numeric columns, just take first 50
        top_genes = df.head(50)

    if len(top_genes) > 0:
        text_parts.append("")
        text_parts.append("Example genes (up to 50):")
        for idx, row in top_genes.iterrows():
            gene_name = str(row.get(gene_column, "")).strip()
            normalized_gene = normalize_gene_identifier(gene_name)
            
            gene_info_parts = [normalized_gene]
            
            if log2fc_col and pd.notna(row.get(log2fc_col)):
                log2fc_val = row.get(log2fc_col)
                gene_info_parts.append(f"log2FC={log2fc_val:.2f}")
            
            if pval_col and pd.notna(row.get(pval_col)):
                pval_val = row.get(pval_col)
                gene_info_parts.append(f"p={pval_val:.4f}")
            
            gene_info = " (".join(gene_info_parts) + (")" if len(gene_info_parts) > 1 else "")
            text_parts.append(f"- {gene_info}")

    return "\n".join(text_parts)


def embed_transcriptomics_dataset(
    page_id: str,
    dataset_name: str,
    genes: Set[str],
    df: pd.DataFrame,
    gene_column: str,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Embed transcriptomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        genes: Set of normalized genes
        df: Full DataFrame with DGE data
        gene_column: Name of the gene column
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    cfg = get_config()

    try:
        # Build text representation
        dataset_text = build_dge_text_representation(
            dataset_name=dataset_name,
            genes=genes,
            df=df,
            gene_column=gene_column,
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][TRANSCRIPTOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_trans_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "omics_type": "Transcriptomics",
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
            "[INGEST][TRANSCRIPTOMICS] Generated %d chunk(s) for transcriptomics dataset %s",
            len(chunks),
            page_id,
        )
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Upserted %d vectors to Pinecone",
            len(vectors),
        )
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Updated Embedding IDs and Last Embedded for dataset %s",
            page_id,
        )

    except Exception as e:
        logger.warning(
            "[INGEST][TRANSCRIPTOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical


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

