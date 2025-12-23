"""Download helpers for LINCS Level 5 data.

This module is intentionally URL/env driven so unit tests can mock HTTP and
production can point to a preferred mirror (CLUE.io, GEO, or a pre-filtered subset).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import requests


def _download(url: str, dest: Path, *, timeout: int = 60) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def download_lincs_level5(output_dir: str, max_compounds: int = 1000) -> str:
    """Download a LINCS Level 5 dataset (GCTX/GCT) and return local path.

    Selection order:
    - If `LINCS_LEVEL5_URL` is set: download that file.
    - Else if `LINCS_LEVEL5_SUBSET_URL` is set: download that (recommended for dev/CI).

    Notes:
    - `max_compounds` is included in the signature to align with ingestion behavior;
      this function may choose a subset URL when max_compounds is small.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    url = os.environ.get("LINCS_LEVEL5_URL")
    subset_url = os.environ.get("LINCS_LEVEL5_SUBSET_URL")

    chosen: Optional[str] = None
    if url:
        chosen = url
    elif subset_url:
        chosen = subset_url

    if not chosen:
        raise RuntimeError(
            "No LINCS download URL configured. Set LINCS_LEVEL5_URL (full) or LINCS_LEVEL5_SUBSET_URL (recommended)."
        )

    filename = os.path.basename(chosen.split("?")[0]) or "lincs_level5.gctx"
    dest = out / filename

    if not dest.exists():
        _download(chosen, dest)

    _ = max_compounds
    return str(dest)


__all__ = ["download_lincs_level5"]


