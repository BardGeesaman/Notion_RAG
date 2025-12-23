"""ClinVar variant_summary downloader + parser."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import requests


DEFAULT_CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


def download_clinvar(dest_dir: str) -> str:
    """Download ClinVar variant_summary.txt (gz) into dest_dir and return path to decompressed txt."""
    out_dir = Path(dest_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    url = os.environ.get("CLINVAR_VARIANT_SUMMARY_URL", DEFAULT_CLINVAR_URL)
    gz_path = out_dir / "variant_summary.txt.gz"
    txt_path = out_dir / "variant_summary.txt"

    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with gz_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    # Decompress
    with gzip.open(gz_path, "rb") as fin, txt_path.open("wb") as fout:
        fout.write(fin.read())

    return str(txt_path)


def _norm_chr(ch: Any) -> Optional[str]:
    if ch is None:
        return None
    s = str(ch).strip()
    if not s:
        return None
    return s.replace("chr", "")


def _norm_allele(a: Any) -> Optional[str]:
    if a is None:
        return None
    s = str(a).strip()
    if not s:
        return None
    return s.upper()


def _as_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    try:
        return int(s)
    except Exception:
        return None


def parse_clinvar(path: str) -> Dict[Tuple[str, int, str, str], Dict[str, Any]]:
    """Parse ClinVar variant_summary.txt into a lookup dict keyed by (chr, pos, ref, alt)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    out: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
    with p.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = _norm_chr(r.get("Chromosome"))
            pos = _as_int(r.get("Start"))
            ref = _norm_allele(r.get("ReferenceAllele"))
            alt = _norm_allele(r.get("AlternateAllele"))
            if chrom is None or pos is None or ref is None or alt is None:
                continue

            var_id = r.get("VariationID") or r.get("VariationID/AlleleID") or r.get("AlleleID")
            clin_sig = r.get("ClinicalSignificance") or r.get("ClinicalSignificanceDescription")
            review_status = r.get("ReviewStatus")
            condition = r.get("PhenotypeList") or r.get("Condition(s)") or r.get("PhenotypeIDS")

            key = (chrom, int(pos), ref, alt)
            out[key] = {
                "clinvar_id": str(var_id).strip() if var_id is not None else None,
                "clinical_significance": str(clin_sig).strip() if clin_sig is not None else None,
                "review_status": str(review_status).strip() if review_status is not None else None,
                "condition": str(condition).strip() if condition is not None else None,
            }

    return out


__all__ = ["download_clinvar", "parse_clinvar"]


