"""AlphaFold DB fetcher."""

from __future__ import annotations

import requests


def fetch_from_alphafold(uniprot_id: str) -> bytes:
    """Fetch an AlphaFold model (v4) from AlphaFold DB as raw bytes (PDB format)."""
    uniprot_id = (uniprot_id or "").strip()
    if not uniprot_id:
        raise ValueError("uniprot_id required")
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.content


__all__ = ["fetch_from_alphafold"]


