"""Parse AutoDock Vina outputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class VinaResult:
    binding_affinity: Optional[float]
    rmsd_lb: Optional[float]
    rmsd_ub: Optional[float]
    pdbqt_content: str


_VINA_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)"
)


def parse_vina_output(pdbqt_path: str, log_path: str) -> VinaResult:
    """Parse output PDBQT and Vina log (log is optional for now)."""
    out_path = Path(pdbqt_path)
    if not out_path.exists():
        raise FileNotFoundError(str(out_path))

    content = out_path.read_text(encoding="utf-8", errors="replace")
    aff = None
    lb = None
    ub = None
    for ln in content.splitlines():
        m = _VINA_RE.search(ln)
        if m:
            aff = float(m.group(1))
            lb = float(m.group(2))
            ub = float(m.group(3))
            break

    _ = log_path  # kept for future parsing; avoids unused argument lint
    return VinaResult(binding_affinity=aff, rmsd_lb=lb, rmsd_ub=ub, pdbqt_content=content)


__all__ = ["VinaResult", "parse_vina_output"]


