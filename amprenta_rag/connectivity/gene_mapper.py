"""NCBI gene metadata mapping utilities (symbol -> Entrez)."""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Dict, List

import requests


NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def _cache_dir() -> Path:
    return Path(os.environ.get("LINCS_CACHE_DIR", "data/lincs"))


def _download(url: str, dest: Path, *, timeout: int = 60) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def load_gene_info(gene_info_path: str) -> Dict[str, int]:
    """Load NCBI gene_info.gz and return mapping {symbol: entrez_id}."""
    p = Path(gene_info_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    out: Dict[str, int] = {}
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            # tax_id, GeneID, Symbol, ...
            try:
                entrez = int(parts[1])
            except Exception:
                continue
            symbol = parts[2].strip()
            if symbol:
                out[symbol] = entrez
    return out


def map_symbols_to_entrez(symbols: List[str]) -> Dict[str, int]:
    """Map gene symbols to Entrez IDs using cached gene_info.gz (download if missing)."""
    cache = _cache_dir()
    gene_info = cache / "Homo_sapiens.gene_info.gz"
    if not gene_info.exists():
        _download(NCBI_GENE_INFO_URL, gene_info)

    mapping = load_gene_info(str(gene_info))
    out: Dict[str, int] = {}
    for s in symbols:
        key = (s or "").strip()
        if key and key in mapping:
            out[key] = mapping[key]
    return out


__all__ = ["load_gene_info", "map_symbols_to_entrez"]


