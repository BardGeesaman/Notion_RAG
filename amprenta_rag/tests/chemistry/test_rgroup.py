from types import SimpleNamespace
from unittest.mock import Mock

import pytest


@pytest.mark.unit
def test_find_common_core_returns_none_when_rdkit_unavailable(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", False)
    assert rgroup.find_common_core(["CCO", "CCN"]) is None


@pytest.mark.unit
def test_find_common_core_empty_list_returns_none():
    from amprenta_rag.chemistry import rgroup

    assert rgroup.find_common_core([]) is None


@pytest.mark.unit
def test_find_common_core_requires_two_valid_molecules(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(rgroup, "_mol_from_smiles_best_effort", lambda smi: None)
    assert rgroup.find_common_core(["BAD1", "BAD2", "BAD3"]) is None


@pytest.mark.unit
def test_find_common_core_returns_none_when_mcs_is_empty(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(rgroup, "_mol_from_smiles_best_effort", lambda smi: {"mol": smi})

    fake_mcs = SimpleNamespace(numAtoms=0, smartsString="")
    monkeypatch.setattr(rgroup, "rdFMCS", SimpleNamespace(FindMCS=Mock(return_value=fake_mcs)))
    assert rgroup.find_common_core(["CCO", "CCN"]) is None


@pytest.mark.unit
def test_find_common_core_success_returns_smarts(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(rgroup, "_mol_from_smiles_best_effort", lambda smi: {"mol": smi})

    fake_mcs = SimpleNamespace(numAtoms=6, smartsString="c1ccccc1")
    monkeypatch.setattr(rgroup, "rdFMCS", SimpleNamespace(FindMCS=Mock(return_value=fake_mcs)))
    assert rgroup.find_common_core(["c1ccccc1", "c1ccncc1"]) == "c1ccccc1"


@pytest.mark.unit
def test_decompose_rgroups_returns_empty_when_rdkit_unavailable(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", False)
    assert rgroup.decompose_rgroups(["CCO"], "core") == []


@pytest.mark.unit
def test_decompose_rgroups_requires_core_smarts(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)
    assert rgroup.decompose_rgroups(["CCO"], "") == []


@pytest.mark.unit
def test_decompose_rgroups_invalid_core_smarts_returns_empty(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)
    fake_chem = SimpleNamespace(MolFromSmarts=Mock(return_value=None))
    monkeypatch.setattr(rgroup, "Chem", fake_chem)
    assert rgroup.decompose_rgroups(["CCO"], "not-a-smarts") == []


@pytest.mark.unit
def test_decompose_rgroups_includes_no_match_error_and_rgroups(monkeypatch):
    from amprenta_rag.chemistry import rgroup

    monkeypatch.setattr(rgroup, "RDKIT_AVAILABLE", True)

    core_pattern = object()

    class FakeAtom:
        def __init__(self, idx: int):
            self._idx = idx
            self._neighbors = []

        def GetIdx(self):
            return self._idx

        def GetNeighbors(self):
            return self._neighbors

    class FakeMol:
        def __init__(self, match, attachment_pairs, has_match=True):
            self._match = tuple(match)
            self._has_match = has_match
            self._atoms = [FakeAtom(i) for i in range(5)]
            for a, b in attachment_pairs:
                self._atoms[a]._neighbors.append(self._atoms[b])
                self._atoms[b]._neighbors.append(self._atoms[a])

        def HasSubstructMatch(self, _pattern):
            return self._has_match

        def GetSubstructMatches(self, _pattern):
            return [self._match] if self._has_match else []

        def GetNumAtoms(self):
            return len(self._atoms)

        def GetAtomWithIdx(self, i: int):
            return self._atoms[i]

    # One molecule matches core with a single attachment point (0 -> 4)
    mol_match = FakeMol(match=(0, 1, 2), attachment_pairs=[(0, 4)], has_match=True)
    # One molecule does not match core
    mol_nomatch = FakeMol(match=(0, 1, 2), attachment_pairs=[], has_match=False)

    def mol_from_smiles(smi: str):
        if smi == "BAD":
            return None
        if smi == "NOMATCH":
            return mol_nomatch
        return mol_match

    fake_chem = SimpleNamespace(
        MolFromSmarts=Mock(return_value=core_pattern),
        PathToSubmol=Mock(return_value={"frag": "x"}),
        MolToSmiles=Mock(return_value="Cl"),
    )
    monkeypatch.setattr(rgroup, "Chem", fake_chem)
    monkeypatch.setattr(rgroup, "_mol_from_smiles_best_effort", mol_from_smiles)

    out = rgroup.decompose_rgroups(["MATCH", "NOMATCH", "BAD"], "core-smarts")
    # BAD is skipped; NOMATCH yields an error row; MATCH yields an R1
    assert out == [
        {"smiles": "MATCH", "R1": "Cl"},
        {"smiles": "NOMATCH", "error": "No match"},
    ]


@pytest.mark.unit
def test_get_rgroup_statistics_empty_returns_empty():
    from amprenta_rag.chemistry import rgroup

    assert rgroup.get_rgroup_statistics([]) == {}


@pytest.mark.unit
def test_get_rgroup_statistics_counts_positions_and_values():
    from amprenta_rag.chemistry import rgroup

    decomposition = [
        {"smiles": "A", "R1": "Cl", "R2": "Br"},
        {"smiles": "B", "R1": "Cl"},
        {"smiles": "C", "R2": "Br"},
        {"smiles": "D", "error": "No match"},
    ]
    stats = rgroup.get_rgroup_statistics(decomposition)
    assert stats == {"R1": {"Cl": 2}, "R2": {"Br": 2}}


