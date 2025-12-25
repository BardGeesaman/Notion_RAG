from __future__ import annotations

import numpy as np
import pytest


pytest.importorskip("rdkit")


def test_generate_conformers_returns_mol_with_conformers():
    from amprenta_rag.chemistry.conformers import generate_conformers

    mols = generate_conformers("CCO", n_conformers=3)
    assert isinstance(mols, list)
    assert len(mols) >= 1
    assert all(m.GetNumConformers() == 1 for m in mols)


def test_optimize_conformer_reduces_energy():
    from amprenta_rag.chemistry.conformers import conformer_energies, generate_conformers, optimize_conformer

    m = generate_conformers("CCO", n_conformers=1)[0]
    e0 = conformer_energies(m, prefer="MMFF")[0]
    optimize_conformer(m, force_field="MMFF")
    e1 = conformer_energies(m, prefer="MMFF")[0]
    # Optimization should not increase energy (allow tiny numerical wiggle).
    assert np.isfinite(e0) and np.isfinite(e1)
    assert e1 <= e0 + 1e-6


def test_get_lowest_energy_conformer():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from amprenta_rag.chemistry.conformers import conformer_energies, get_lowest_energy_conformer

    m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    try:
        params.maxAttempts = 50  # type: ignore[attr-defined]
    except Exception:
        pass
    conf_ids = list(AllChem.EmbedMultipleConfs(m, numConfs=5, params=params))
    assert conf_ids

    # Optimize so energies are meaningful
    for cid in conf_ids:
        try:
            AllChem.MMFFOptimizeMolecule(m, confId=int(cid))
        except Exception:
            pass

    energies = conformer_energies(m, prefer="MMFF")
    min_e = float(np.nanmin(np.asarray(energies, dtype=float)))

    best = get_lowest_energy_conformer(m)
    best_e = conformer_energies(best, prefer="MMFF")[0]
    assert best.GetNumConformers() == 1
    assert abs(best_e - min_e) < 1e-3


def test_conformer_to_pdb_returns_string():
    from amprenta_rag.chemistry.conformers import conformer_to_pdb, generate_conformers

    m = generate_conformers("CCO", n_conformers=1)[0]
    pdb = conformer_to_pdb(m, conf_id=0)
    assert isinstance(pdb, str)
    assert "ATOM" in pdb or "HETATM" in pdb


def test_align_conformers_no_error():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from amprenta_rag.chemistry.conformers import align_conformers

    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    try:
        params.maxAttempts = 50  # type: ignore[attr-defined]
    except Exception:
        pass
    conf_ids = list(AllChem.EmbedMultipleConfs(m, numConfs=3, params=params))
    assert len(conf_ids) >= 2
    # Should not raise
    align_conformers(m)


def test_invalid_smiles_raises_error():
    from amprenta_rag.chemistry.conformers import generate_conformers

    with pytest.raises(ValueError):
        generate_conformers("not-a-smiles", n_conformers=1)


