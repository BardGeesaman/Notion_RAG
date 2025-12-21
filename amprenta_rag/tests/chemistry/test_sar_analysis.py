from __future__ import annotations

import sys
from types import SimpleNamespace, ModuleType
from uuid import UUID

import pandas as pd

from amprenta_rag.chemistry import sar_analysis as sa


class _FakeDBGen:
    def __init__(self, session):
        self.session = session
        self.used = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.used:
            raise StopIteration
        self.used = True
        return self.session

    def close(self):
        if hasattr(self.session, "close"):
            self.session.close()


class _FakeQuery:
    def __init__(self, rows, allowed_ids=None):
        self.rows = rows
        self._limit = None
        self.allowed_ids = allowed_ids

    def limit(self, n):
        self._limit = n
        return self

    def filter(self, *args, **kwargs):
        if self.allowed_ids is not None:
            self.rows = [r for r in self.rows if r.compound_id in self.allowed_ids]
        return self

    def all(self):
        if self._limit is not None:
            return self.rows[: self._limit]
        return self.rows


class _FakeSession:
    def __init__(self, rows, allowed_ids=None):
        self.rows = rows
        self.allowed_ids = allowed_ids
        self.closed = False

    def query(self, model):
        return _FakeQuery(self.rows, allowed_ids=self.allowed_ids)

    def close(self):
        self.closed = True


class _FakeCompound:
    def __init__(self, cid: str, smiles: str, mw: float = 300.0):
        self.compound_id = cid
        self.smiles = smiles
        self.molecular_weight = mw
        self.logp = 2.0
        self.hbd_count = 1
        self.hba_count = 2
        self.rotatable_bonds = 3


def test_get_compound_properties(monkeypatch):
    rows = [_FakeCompound("c1", "C"), _FakeCompound("c2", "CC")]
    sess = _FakeSession(rows)
    monkeypatch.setattr(sa, "get_db", lambda: _FakeDBGen(sess))
    df = sa.get_compound_properties(limit=10)
    assert set(df["compound_id"]) == {"c1", "c2"}

    # Empty path
    sess_empty = _FakeSession([])
    monkeypatch.setattr(sa, "get_db", lambda: _FakeDBGen(sess_empty))
    assert sa.get_compound_properties().empty


def test_get_activity_data_filters(monkeypatch):
    rows = [_FakeCompound("c1", "C"), _FakeCompound("c2", "CC")]
    sess = _FakeSession(rows, allowed_ids=["c2"])
    monkeypatch.setattr(sa, "get_db", lambda: _FakeDBGen(sess))
    df = sa.get_activity_data(compound_ids=["c2"])
    assert list(df["compound_id"]) == ["c2"]


def test_calculate_lipinski_with_fake_rdkit(monkeypatch):
    fake_rdkit = ModuleType("rdkit")

    class FakeMol:
        def __init__(self, smiles):
            self.smiles = smiles

    class FakeChem(ModuleType):
        @staticmethod
        def MolFromSmiles(smiles):
            return FakeMol(smiles)

    class FakeDescriptors(ModuleType):
        @staticmethod
        def MolWt(mol):
            return 100.0

        @staticmethod
        def MolLogP(mol):
            return 1.0

        @staticmethod
        def NumHDonors(mol):
            return 1

        @staticmethod
        def NumHAcceptors(mol):
            return 2

    fake_rdkit.Chem = FakeChem("Chem")
    fake_rdkit.Chem.Descriptors = FakeDescriptors("Descriptors")
    fake_descriptors = FakeDescriptors("Descriptors")

    monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_rdkit.Chem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.Descriptors", fake_descriptors)

    res = sa.calculate_lipinski("C")
    assert res["valid"] is True
    assert res["passes_ro5"] is True


def test_detect_activity_cliffs_with_fake_rdkit(monkeypatch):
    # Build fake rdkit modules
    fake_rdkit = ModuleType("rdkit")

    class FakeMol:
        def __init__(self, smiles):
            self.smiles = smiles

    class FakeChem(ModuleType):
        @staticmethod
        def MolFromSmiles(smiles):
            return FakeMol(smiles)

    class FakeAllChem(ModuleType):
        @staticmethod
        def GetMorganFingerprintAsBitVect(mol, radius, nBits=2048):
            return f"fp-{mol.smiles}"

    class FakeDataStructs(ModuleType):
        @staticmethod
        def TanimotoSimilarity(a, b):
            return 0.8

    fake_rdkit.Chem = FakeChem("Chem")
    fake_rdkit.Chem.AllChem = FakeAllChem("AllChem")
    fake_rdkit.Chem.DataStructs = FakeDataStructs("DataStructs")
    fake_rdkit.DataStructs = fake_rdkit.Chem.DataStructs
    monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_rdkit.Chem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.AllChem", fake_rdkit.Chem.AllChem)
    monkeypatch.setitem(sys.modules, "rdkit.Chem.DataStructs", fake_rdkit.Chem.DataStructs)
    monkeypatch.setitem(sys.modules, "rdkit.DataStructs", fake_rdkit.DataStructs)

    # Fake ActivityResult and Compound rows
    class FakeActivity:
        def __init__(self, cid, smiles, assay_id, value):
            self.compound = _FakeCompound(cid, smiles)
            self.compound_id = uuid4()
            self.assay_id = UUID(assay_id)
            self.value = value

    results = [
        FakeActivity("c1", "C", assay_id="00000000-0000-0000-0000-000000000001", value=1.0),
        FakeActivity("c2", "CC", assay_id="00000000-0000-0000-0000-000000000001", value=20.0),
    ]

    class FakeQuery:
        def __init__(self, rows):
            self.rows = rows

        def join(self, *_):
            return self

        def filter(self, *_):
            return self

        def all(self):
            return self.rows

    class FakeSession:
        def __init__(self, rows):
            self.rows = rows

        def query(self, model):
            return FakeQuery(self.rows)

        def close(self):
            pass

    monkeypatch.setattr(sa, "get_db", lambda: _FakeDBGen(FakeSession(results)))

    cliffs = sa.detect_activity_cliffs()
    assert cliffs
    assert cliffs[0]["fold_change"] >= 10.0

