"""Ligand preparation utilities (SMILES -> 3D -> PDBQT)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path


def prepare_ligand(smiles: str, output_pdbqt: str) -> bool:
    """Prepare a ligand from SMILES and write a PDBQT file.

    Steps:
    - RDKit: parse SMILES, add Hs, embed 3D (ETKDG), MMFF optimize
    - Convert to PDBQT:
      - Prefer Meeko if installed (pure-Python)
      - Fallback to OpenBabel CLI (`obabel`) if available
    """
    out_path = Path(output_pdbqt)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return False

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, params) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        # Some molecules may not have MMFF params; still proceed
        pass

    # Prefer Meeko if available
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            return False
        pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            _ = err
            return False
        out_path.write_text(pdbqt_str, encoding="utf-8")
        return True
    except Exception:
        pass

    # Fallback: OpenBabel CLI
    if shutil.which("obabel") is None:
        return False

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        pdb_path = td_path / "ligand.pdb"
        Chem.MolToPDBFile(mol, str(pdb_path))
        cmd = ["obabel", "-i", "pdb", str(pdb_path), "-o", "pdbqt", "-O", str(out_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return proc.returncode == 0 and out_path.exists()


__all__ = ["prepare_ligand"]


