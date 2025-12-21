from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.utils import batch_ops


class FakeDataset:
    def __init__(self, name: str):
        self.id = uuid4()
        self.name = name
        self.omics_type = "rna"
        self.description = None
        self.organism = ["human"]
        self.disease = ["flu"]
        self.created_at = None


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)
        self._delete_count = len(self._data)

    def filter(self, *args, **kwargs):
        return self

    def all(self):
        return list(self._data)

    def delete(self, synchronize_session=False):
        return self._delete_count


class FakeSession:
    def __init__(self, datasets=None):
        self.datasets = datasets or []
        self.deleted = 0
        self.committed = False

    def query(self, model):
        if model is batch_ops.Dataset:
            return FakeQuery(self.datasets)
        return FakeQuery([])

    def commit(self):
        self.committed = True


def test_batch_export_datasets_json(monkeypatch):
    db = FakeSession([FakeDataset("ds1")])
    monkeypatch.setattr(batch_ops, "export_experiments", lambda ids, fmt, db: b"exp")
    monkeypatch.setattr(batch_ops, "export_compounds", lambda ids, fmt, db: b"cmp")
    monkeypatch.setattr(batch_ops, "export_signatures", lambda ids, fmt, db: b"sig")
    monkeypatch.setattr(batch_ops, "export_to_csv", lambda df: b"csv", raising=False)
    monkeypatch.setattr(batch_ops, "export_to_json", lambda df: b"json", raising=False)
    monkeypatch.setattr(batch_ops, "export_to_excel", lambda df: b"xlsx", raising=False)

    out = batch_ops.batch_export("dataset", [str(db.datasets[0].id)], db, format="json")
    assert b"ds1" in out


def test_batch_export_unsupported():
    with pytest.raises(ValueError):
        batch_ops.batch_export("unknown", [], None)


def test_batch_delete_counts(monkeypatch):
    db = FakeSession()
    db.deleted = 3

    def fake_delete(*args, **kwargs):
        return db.deleted

    class E:
        id = type("col", (), {"in_": staticmethod(lambda ids: ("in", ids))})

    class C(E):
        pass

    class S(E):
        pass

    class D(E):
        pass

    monkeypatch.setattr(batch_ops, "Experiment", E)
    monkeypatch.setattr(batch_ops, "Compound", C)
    monkeypatch.setattr(batch_ops, "Signature", S)
    monkeypatch.setattr(batch_ops, "Dataset", D)

    class Q:
        def __init__(self, deleted):
            self.deleted = deleted

        def filter(self, *a, **k):
            return self

        def delete(self, synchronize_session=False):
            return self.deleted

    class Session:
        def query(self, model):
            return Q(db.deleted)

        def commit(self):
            db.committed = True

    session = Session()
    deleted = batch_ops.batch_delete("experiment", [str(uuid4())], session)
    assert deleted == db.deleted
    assert db.committed is True

