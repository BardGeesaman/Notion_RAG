"""VEP (--tab) TSV parser."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional


def _parse_location(loc: str) -> tuple[Optional[str], Optional[int]]:
    """Parse VEP Location field like '1:12345' or '1:12345-12345'."""
    if not loc:
        return None, None
    t = str(loc).strip()
    if ":" not in t:
        return None, None
    chrom, rest = t.split(":", 1)
    chrom = chrom.replace("chr", "").strip()
    if "-" in rest:
        rest = rest.split("-", 1)[0]
    try:
        pos = int(rest)
    except Exception:
        pos = None
    return chrom or None, pos


def _as_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s in (".", "-", "NA", "nan"):
        return None
    try:
        return float(s)
    except Exception:
        return None


def _parse_pred(x: Any) -> Optional[str]:
    """Parse strings like 'deleterious(0.02)' into 'deleterious'."""
    if x is None:
        return None
    s = str(x).strip()
    if not s or s in (".", "-", "NA"):
        return None
    if "(" in s:
        return s.split("(", 1)[0].strip() or None
    return s


def _get(row: Dict[str, Any], *names: str) -> Any:
    for n in names:
        if n in row:
            return row.get(n)
    return None


def parse_vep_tsv(path: str) -> List[Dict[str, Any]]:
    """Parse VEP --tab output TSV into a list of dicts.

    Extracts:
    - Location (chr:pos)
    - Allele
    - Gene
    - Consequence
    - IMPACT
    - SIFT
    - PolyPhen
    - gnomAD_AF
    - CADD_PHRED
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    rows: List[Dict[str, Any]] = []
    with p.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            loc = _get(r, "Location", "LOCATION")
            allele = _get(r, "Allele", "ALLELE")
            gene = _get(r, "Gene", "GENE")
            consequence = _get(r, "Consequence", "CONSEQUENCE")
            impact = _get(r, "IMPACT", "Impact")
            sift = _get(r, "SIFT", "SIFT_prediction")
            polyphen = _get(r, "PolyPhen", "PolyPhen_prediction")
            gnomad_af = _get(r, "gnomAD_AF", "gnomAD_AF_1KG", "gnomADg_AF", "gnomad_af")
            cadd = _get(r, "CADD_PHRED", "CADD_PHRED_score", "CADD_PHRED_PHRED", "CADD_PHRED")

            chrom, pos = _parse_location(str(loc) if loc is not None else "")

            # Optional symbol field if present
            symbol = _get(r, "SYMBOL", "Symbol", "HGNC", "HGNC_ID")

            rows.append(
                {
                    "chromosome": chrom,
                    "position": pos,
                    "alt_allele": str(allele).strip() if allele is not None else None,
                    "gene": str(gene).strip() if gene is not None else None,
                    "gene_symbol": str(symbol).strip() if symbol is not None else None,
                    "consequence": str(consequence).strip() if consequence is not None else None,
                    "impact": str(impact).strip() if impact is not None else None,
                    "sift_prediction": _parse_pred(sift),
                    "polyphen_prediction": _parse_pred(polyphen),
                    "gnomad_af": _as_float(gnomad_af),
                    "cadd_phred": _as_float(cadd),
                    "raw": r,
                }
            )
    return rows


__all__ = ["parse_vep_tsv"]


