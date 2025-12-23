from __future__ import annotations

import pytest


def test_validate_smiles():
    pytest.importorskip("rdkit")
    from amprenta_rag.analysis.sar_whatif import validate_smiles

    assert validate_smiles("c1ccccc1") is True
    assert validate_smiles("not_a_smiles") is False


def test_scaffold_hop_benzene_to_pyridine():
    pytest.importorskip("rdkit")
    from amprenta_rag.analysis.sar_whatif import scaffold_hop

    out = scaffold_hop("c1ccccc1", "benzene_to_pyridine")
    assert out
    assert any(("n" in s) and ("c1" in s) for s in out)


