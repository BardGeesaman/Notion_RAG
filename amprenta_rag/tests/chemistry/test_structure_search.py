from __future__ import annotations

import types

from amprenta_rag.chemistry import structure_search as ss


def test_substructure_search_returns_empty_when_no_rdkit(monkeypatch):
    monkeypatch.setattr(ss, "RDKIT_AVAILABLE", False)
    assert ss.substructure_search("C1=CC=CC=C1") == []


def test_similarity_search_returns_empty_when_no_rdkit(monkeypatch):
    monkeypatch.setattr(ss, "RDKIT_AVAILABLE", False)
    assert ss.similarity_search("C", threshold=0.5) == []


def test_substructure_search_invalid_smarts(monkeypatch):
    # Force RDKit to be "available" with simple parser
    fake_chem = types.SimpleNamespace(
        MolFromSmarts=lambda s: None,
        MolFromSmiles=lambda s: None,
    )
    monkeypatch.setattr(ss, "Chem", fake_chem)
    monkeypatch.setattr(ss, "RDKIT_AVAILABLE", True)
    # Should return [] due to pattern None
    assert ss.substructure_search("not-a-smarts") == []

