"""Parse fpocket output directories."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass(frozen=True)
class PocketInfo:
    pocket_rank: int
    score: Optional[float] = None
    volume: Optional[float] = None
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    center_z: Optional[float] = None
    residues: Optional[List[str]] = None
    pocket_pdb_path: Optional[str] = None


_FLOAT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)")
_REMARK_RES_TRIPLE = re.compile(r"([A-Za-z0-9]):([A-Za-z]{3}):(-?\d+)")
_REMARK_RES_SPACED = re.compile(r"\b([A-Za-z0-9])\s+([A-Za-z]{3})\s+(-?\d+)\b")


def _parse_residues_from_pocket_pdb(pdb_path: Path) -> Optional[List[str]]:
    """Best-effort extraction of residues from pocket PDB REMARK lines.

    fpocket output varies by build/config; for MVP we only parse explicit residue
    annotations in REMARK lines (if present). Otherwise returns None.
    """
    residues: set[str] = set()
    try:
        for ln in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not ln.startswith("REMARK"):
                continue
            for m in _REMARK_RES_TRIPLE.finditer(ln):
                residues.add(f"{m.group(1)}:{m.group(2).upper()}:{int(m.group(3))}")
            # Some variants may store "A GLY 12" patterns in REMARKs.
            for m in _REMARK_RES_SPACED.finditer(ln):
                residues.add(f"{m.group(1)}:{m.group(2).upper()}:{int(m.group(3))}")
    except Exception:
        return None

    return sorted(residues) if residues else None


def _parse_info_txt(path: Path) -> List[PocketInfo]:
    # Supports typical fpocket *_info.txt formatting, but stays permissive.
    pockets: List[PocketInfo] = []
    cur_rank: Optional[int] = None
    cur_score: Optional[float] = None
    cur_vol: Optional[float] = None

    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Pocket" in ln and "Pocket" == ln.strip().split()[0]:
            # Pocket 1
            m = re.search(r"Pocket\s+(\d+)", ln)
            if m:
                # Flush previous
                if cur_rank is not None:
                    pockets.append(PocketInfo(pocket_rank=cur_rank, score=cur_score, volume=cur_vol))
                cur_rank = int(m.group(1))
                cur_score = None
                cur_vol = None
            continue
        if "Score" in ln and ":" in ln:
            m = _FLOAT_RE.search(ln)
            if m:
                cur_score = float(m.group(1))
        if "Volume" in ln and ":" in ln:
            m = _FLOAT_RE.search(ln)
            if m:
                cur_vol = float(m.group(1))

    if cur_rank is not None:
        pockets.append(PocketInfo(pocket_rank=cur_rank, score=cur_score, volume=cur_vol))
    return pockets


def _parse_centers_from_pqr(path: Path) -> dict[int, tuple[float, float, float]]:
    # Heuristic: compute centroid for each pocket based on ATOM coordinates.
    # We assume blocks separated by "TER" correspond to pockets in order.
    centers: dict[int, tuple[float, float, float]] = {}
    coords: List[tuple[float, float, float]] = []
    pocket_idx = 1

    def _flush():
        nonlocal coords, pocket_idx
        if coords:
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            zs = [c[2] for c in coords]
            centers[pocket_idx] = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
            pocket_idx += 1
            coords = []

    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith("TER"):
            _flush()
            continue
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            parts = ln.split()
            # Try last 3 floats
            floats = []
            for p in parts:
                try:
                    floats.append(float(p))
                except Exception:
                    continue
            if len(floats) >= 3:
                coords.append((floats[-3], floats[-2], floats[-1]))
    _flush()
    return centers


def parse_fpocket_output(output_dir: str) -> List[PocketInfo]:
    """Parse fpocket output and return top 5 pockets sorted by score."""
    out = Path(output_dir)
    if not out.exists():
        raise FileNotFoundError(str(out))

    info_files = list(out.glob("*_info.txt")) + list(out.glob("*.info")) + list(out.glob("pockets_info.txt"))
    if not info_files:
        # try inside fpocket output folder
        info_files = list(out.rglob("*_info.txt"))
    if not info_files:
        return []

    pockets = _parse_info_txt(info_files[0])

    pqr_files = list(out.glob("*_pockets.pqr")) + list(out.rglob("*_pockets.pqr"))
    centers = _parse_centers_from_pqr(pqr_files[0]) if pqr_files else {}

    # Attach centers + pocket PDB path if present
    for i, p in enumerate(pockets):
        c = centers.get(p.pocket_rank)
        pocket_pdb = None
        residues = None
        # Common: pocket{rank}.pdb under *_out/pockets/
        cand = list(out.rglob(f"pocket{p.pocket_rank}.pdb"))
        if cand:
            pocket_pdb = str(cand[0])
            residues = _parse_residues_from_pocket_pdb(Path(pocket_pdb))
        pockets[i] = PocketInfo(
            pocket_rank=p.pocket_rank,
            score=p.score,
            volume=p.volume,
            center_x=c[0] if c else None,
            center_y=c[1] if c else None,
            center_z=c[2] if c else None,
            residues=residues,
            pocket_pdb_path=pocket_pdb,
        )

    pockets.sort(key=lambda x: float(x.score or 0.0), reverse=True)
    return pockets[:5]


__all__ = ["PocketInfo", "parse_fpocket_output"]


