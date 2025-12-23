"""Protein structure preparation utilities (PDBFixer + OpenMM).

This module is intentionally optional-dependency based: environments without
`pdbfixer`/`openmm` can still import the package, but preparation will raise
a helpful ImportError.
"""

from __future__ import annotations

from pathlib import Path


def prepare_structure(input_pdb: str, output_pdb: str, chain_id: str | None = None) -> None:
    """Prepare a PDB file for downstream docking/simulation.

    Steps:
    - Optional: keep only a single chain (if chain_id provided)
    - findMissingResidues
    - findMissingAtoms
    - addMissingAtoms
    - addMissingHydrogens(pH=7.4)
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "Structure preparation requires `pdbfixer` and `openmm`. "
            "Install with: pip install pdbfixer openmm"
        ) from e

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


