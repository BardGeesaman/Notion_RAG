"""Protein structure preparation utilities (PDBFixer + OpenMM).

This module is intentionally optional-dependency based: environments without
`pdbfixer`/`openmm` can still import the package, but preparation will raise
a helpful ImportError.
"""

from __future__ import annotations

from pathlib import Path

try:
    from pdbfixer import PDBFixer  # noqa: F401
    import openmm.app  # noqa: F401
    STRUCTURAL_AVAILABLE = True
except ImportError:
    STRUCTURAL_AVAILABLE = False


def _check_structural_deps():
    """Check if structural biology dependencies are available."""
    if not STRUCTURAL_AVAILABLE:
        raise ImportError(
            "Structural biology dependencies not installed. "
            "Install with: pip install -r requirements-structural.txt "
            "or: conda install -c conda-forge pdbfixer openmm"
        )


def prepare_structure(input_pdb: str, output_pdb: str, chain_id: str | None = None) -> None:
    """Prepare a PDB file for downstream docking/simulation.

    Steps:
    - Optional: keep only a single chain (if chain_id provided)
    - findMissingResidues
    - findMissingAtoms
    - addMissingAtoms
    - addMissingHydrogens(pH=7.4)
    """
    _check_structural_deps()
    
    # Import after checking dependencies
    from openmm.app import PDBFile

    in_path = Path(input_pdb)
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))

    out_path = Path(output_pdb)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fixer = PDBFixer(filename=str(in_path))

    if chain_id:
        cid = str(chain_id).strip()
        if cid:
            keep = {cid}
            remove = [c.id for c in fixer.topology.chains() if c.id not in keep]
            if remove:
                fixer.removeChains(remove)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    with out_path.open("w", encoding="utf-8") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)


__all__ = ["prepare_structure"]


