"""MAGeCK result parsing utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import pandas as pd


def parse_gene_summary(path: str) -> List[Dict[str, Any]]:
    """Parse MAGeCK *.gene_summary.txt into a list of dicts.

    Returned keys:
    - gene
    - neg_lfc
    - pos_lfc
    - neg_p
    - pos_p
    - fdr
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    df = pd.read_csv(p, sep="\t")
    cols = {c.strip(): c for c in df.columns}

    # MAGeCK uses "id" for gene symbol (common) and column names like "neg|lfc"
    gene_col = cols.get("id") or cols.get("gene") or cols.get("Gene") or cols.get("symbol")
    if not gene_col:
        raise ValueError("gene_summary missing gene identifier column (expected id/gene/Gene/symbol)")

    def _get(*names: str) -> str | None:
        for n in names:
            if n in cols:
                return cols[n]
        return None

    neg_lfc_col = _get("neg|lfc", "neg_lfc", "neg.lfc")
    pos_lfc_col = _get("pos|lfc", "pos_lfc", "pos.lfc")
    neg_p_col = _get("neg|p-value", "neg_p", "neg.p", "neg_p_value")
    pos_p_col = _get("pos|p-value", "pos_p", "pos.p", "pos_p_value")
    fdr_col = _get("fdr", "FDR")

    out: List[Dict[str, Any]] = []
    for _, r in df.iterrows():
        out.append(
            {
                "gene": r.get(gene_col),
                "neg_lfc": float(r.get(neg_lfc_col)) if neg_lfc_col and pd.notna(r.get(neg_lfc_col)) else None,
                "pos_lfc": float(r.get(pos_lfc_col)) if pos_lfc_col and pd.notna(r.get(pos_lfc_col)) else None,
                "neg_p": float(r.get(neg_p_col)) if neg_p_col and pd.notna(r.get(neg_p_col)) else None,
                "pos_p": float(r.get(pos_p_col)) if pos_p_col and pd.notna(r.get(pos_p_col)) else None,
                "fdr": float(r.get(fdr_col)) if fdr_col and pd.notna(r.get(fdr_col)) else None,
            }
        )
    return out


__all__ = ["parse_gene_summary"]


