"""
Text representation building for transcriptomics datasets.

This module provides functions for building text representations of DGE datasets
for RAG embedding, including summary statistics and top genes.
"""

from __future__ import annotations

from typing import List, Optional, Set

import pandas as pd

from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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

