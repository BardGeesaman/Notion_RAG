"""Spectral library parser (MGF)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple


def parse_mgf(path: str) -> List[Dict[str, Any]]:
    """Parse an MGF file into a list of spectra dicts.

    Returns items like:
    {
      "precursor_mz": float,
      "name": str,
      "precursor_type": str|None,
      "collision_energy": float|None,
      "peaks": [(mz, intensity), ...],
    }
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    out: List[Dict[str, Any]] = []
    in_block = False
    meta: Dict[str, str] = {}
    peaks: List[Tuple[float, float]] = []

    def _flush():
        nonlocal meta, peaks
        if not meta and not peaks:
            return
        precursor_mz = None
        pepmass = meta.get("PEPMASS")
        if pepmass:
            # PEPMASS can be "mz" or "mz intensity"
            try:
                precursor_mz = float(pepmass.split()[0])
            except Exception:
                precursor_mz = None
        name = meta.get("NAME") or meta.get("TITLE") or "unknown"
        precursor_type = meta.get("CHARGE") or meta.get("PRECURSORTYPE")
        ce = None
        if meta.get("COLLISIONENERGY"):
            try:
                ce = float(str(meta.get("COLLISIONENERGY")).split()[0])
            except Exception:
                ce = None
        if precursor_mz is not None and peaks:
            out.append(
                {
                    "precursor_mz": float(precursor_mz),
                    "name": str(name),
                    "precursor_type": str(precursor_type) if precursor_type else None,
                    "collision_energy": ce,
                    "peaks": peaks,
                }
            )
        meta = {}
        peaks = []

    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        t = ln.strip()
        if not t:
            continue
        if t.upper() == "BEGIN IONS":
            in_block = True
            meta = {}
            peaks = []
            continue
        if t.upper() == "END IONS":
            _flush()
            in_block = False
            continue
        if not in_block:
            continue
        if "=" in t:
            k, v = t.split("=", 1)
            meta[k.strip().upper()] = v.strip()
            continue
        # Peak line: mz intensity
        parts = t.split()
        if len(parts) >= 2:
            try:
                mz = float(parts[0])
                inten = float(parts[1])
                peaks.append((mz, inten))
            except Exception:
                continue

    # flush if file ended unexpectedly
    if in_block:
        _flush()

    return out


__all__ = ["parse_mgf"]


