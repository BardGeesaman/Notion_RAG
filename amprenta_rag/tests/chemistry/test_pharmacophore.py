from __future__ import annotations

from amprenta_rag.chemistry import pharmacophore


def test_get_pharmacophore_features_rdkit_missing(monkeypatch):
    monkeypatch.setattr(pharmacophore, "RDKIT_AVAILABLE", False)
    feats = pharmacophore.get_pharmacophore_features("CCO")
    assert feats == {"hbd": 0, "hba": 0, "aromatic_rings": 0, "hydrophobic": 0.0}


def test_pharmacophore_search_rdkit_missing(monkeypatch):
    monkeypatch.setattr(pharmacophore, "RDKIT_AVAILABLE", False)
    assert pharmacophore.pharmacophore_search(min_hbd=1) == []

