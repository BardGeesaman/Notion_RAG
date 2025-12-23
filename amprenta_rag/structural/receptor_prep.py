"""Receptor preparation utilities (ProteinStructure -> PDBQT)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from amprenta_rag.database.models import ProteinStructure


def prepare_receptor(structure: ProteinStructure, output_pdbqt: str) -> bool:
    """Prepare a receptor PDBQT from a ProteinStructure.

    Uses the best available structure file:
    - prefer StructureFile.file_type == "prepared"
    - else fall back to "pdb"

    Conversion is performed via OpenBabel CLI (`obabel`) if available.
    """
    out_path = Path(output_pdbqt)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    files = list(structure.files or [])
    src = next((f for f in files if f.file_type == "prepared"), None) or next((f for f in files if f.file_type == "pdb"), None)
    if src is None:
        return False

    in_path = Path(src.file_path)
    if not in_path.exists():
        return False

    if shutil.which("obabel") is None:
        return False

    # OpenBabel generally expects explicit input format; we use pdb -> pdbqt.
    # Some builds may require flags for receptor (-xr); we keep it minimal for MVP.
    cmd = ["obabel", "-i", "pdb", str(in_path), "-o", "pdbqt", "-O", str(out_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=str(out_path.parent))
    return proc.returncode == 0 and out_path.exists()


__all__ = ["prepare_receptor"]


