"""RCSB PDB fetcher."""

from __future__ import annotations

import requests


def fetch_from_pdb(pdb_id: str) -> bytes:
    """Fetch a structure from the RCSB PDB as raw bytes (PDB format)."""
    pdb_id = (pdb_id or "").strip()
    if not pdb_id:
        raise ValueError("pdb_id required")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.content


__all__ = ["fetch_from_pdb"]


