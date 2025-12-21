from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace, ModuleType
import sys
from uuid import UUID, uuid4

import pytest

from amprenta_rag.ingestion import postgres_signature_creation as psc


class _Expr:
    def __or__(self, other: object):
        return self


class _Col:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):
        return _Expr()


class FakeSignatureModel:
    # Class-level "columns" for filter expressions
    name = _Col("name")
    short_id = _Col("short_id")
    id = _Col("id")

    def __init__(
        self,
        name: str,
        description: str | None = None,
        modalities: list[str] | None = None,
        short_id: str | None = None,
        biomarker_role: list[str] | None = None,
        phenotype_axes: list[str] | None = None,
        data_ownership: str | None = None,
        external_ids: dict | None = None,
    ):
        self.name = name
        self.description = description
        self.modalities = modalities or []
        self.short_id = short_id
        self.biomarker_role = biomarker_role or []
        self.phenotype_axes = phenotype_axes or []
        self.data_ownership = data_ownership
        self.external_ids = external_ids
        self.id = uuid4()


@dataclass
class FakeFeature:
    id: UUID = field(default_factory=uuid4)


class FakeComponentModel:
    signature_id = _Col("signature_id")

    def __init__(
        self,
        signature_id: UUID,
        feature_id: UUID | None,
        feature_name: str,
        feature_type: str,
        direction: str | None,
        weight: float,
    ):
        self.signature_id = signature_id
        self.feature_id = feature_id
        self.feature_name = feature_name
        self.feature_type = feature_type
        self.direction = direction
        self.weight = weight


@dataclass
class FakeExistingComponent:
    feature_name: str
    feature_type: str
    direction: str | None = None


class FakeQuery:
    def __init__(self, first=None, all_=None):
        self._first = first
        self._all = list(all_ or [])

    def filter(self, *_a, **_k):
        return self

    def first(self):
        return self._first

    def all(self):
        return list(self._all)


class FakeDB:
    def __init__(self):
        self.added = []
        self.commits = 0
        self.refreshed = []
        self.rolled_back = 0
        self._queries = {}

    def set_query(self, model, *, first=None, all_=None):
        self._queries[model] = FakeQuery(first=first, all_=all_)

    def query(self, model):
        return self._queries.get(model, FakeQuery(first=None, all_=[]))

    def add(self, obj):
        self.added.append(obj)

    def commit(self):
        self.commits += 1

    def refresh(self, obj):
        self.refreshed.append(obj)

    def rollback(self):
        self.rolled_back += 1


def test_create_signature_new(monkeypatch):
    db = FakeDB()
    monkeypatch.setattr(psc, "SignatureModel", FakeSignatureModel)
    monkeypatch.setattr(psc, "generate_signature_short_id", lambda name, version: "SID")
    db.set_query(FakeSignatureModel, first=None)

    sig = SimpleNamespace(name="SigA", description="D", modalities=["lipid"], components=[])
    created = psc.create_signature_in_postgres(sig, description="DESC", biomarker_roles=["r"], phenotype_axes=["a"], db=db)

    assert isinstance(created, FakeSignatureModel)
    assert created.name == "SigA"
    assert created.short_id == "SID"
    assert created.description == "DESC"
    assert created.biomarker_role == ["r"]
    assert created.phenotype_axes == ["a"]
    assert db.commits == 1
    assert db.refreshed == [created]


def test_create_signature_existing_updates_metadata(monkeypatch):
    existing = FakeSignatureModel(name="SigA", description=None, modalities=[], short_id=None, data_ownership=None)
    db = FakeDB()
    monkeypatch.setattr(psc, "SignatureModel", FakeSignatureModel)
    monkeypatch.setattr(psc, "generate_signature_short_id", lambda name, version: "SID")
    db.set_query(FakeSignatureModel, first=existing)

    sig = SimpleNamespace(name="SigA", description="D", modalities=["gene"], components=[])
    got = psc.create_signature_in_postgres(
        sig,
        description="DESC",
        biomarker_roles=["r"],
        phenotype_axes=["a"],
        data_ownership="Public",
        db=db,
    )

    assert got is existing
    assert existing.description == "DESC"
    assert existing.biomarker_role == ["r"]
    assert existing.phenotype_axes == ["a"]
    assert existing.data_ownership == "Public"
    assert existing.modalities == ["gene"]
    assert existing.short_id == "SID"
    assert db.commits == 1


def test_create_signature_components_creates_and_skips_existing(monkeypatch):
    db = FakeDB()
    sig_model = FakeSignatureModel(name="SigA")
    monkeypatch.setattr(psc, "SignatureComponentModel", FakeComponentModel)
    db.set_query(FakeComponentModel, all_=[FakeExistingComponent(feature_name="A", feature_type="gene", direction=None)])
    monkeypatch.setattr(psc, "find_or_create_feature_in_postgres", lambda name, ftype, db: FakeFeature())

    # One existing (should skip), one new, one with missing name (skip)
    comp_existing = SimpleNamespace(feature_type="gene", feature_name="A", direction=None, weight=1.0)
    comp_new = SimpleNamespace(feature_type="gene", feature_name="B", direction=None, weight=2.0)
    comp_blank = SimpleNamespace(feature_type="gene", feature_name="", direction=None, weight=1.0)
    sig = SimpleNamespace(components=[comp_existing, comp_new, comp_blank])

    count, created = psc.create_signature_components_in_postgres(sig_model, sig, db=db)
    assert count == 1
    assert ("gene", "B") in created
    assert db.commits == 1
    assert len(db.added) == 1
    assert isinstance(db.added[0], FakeComponentModel)


def test_create_signature_components_per_component_error_continues(monkeypatch):
    db = FakeDB()
    sig_model = FakeSignatureModel(name="SigA")
    monkeypatch.setattr(psc, "SignatureComponentModel", FakeComponentModel)
    db.set_query(FakeComponentModel, all_=[])

    def fake_find(name, ftype, db):
        if name == "BAD":
            raise RuntimeError("boom")
        return FakeFeature()

    monkeypatch.setattr(psc, "find_or_create_feature_in_postgres", fake_find)

    sig = SimpleNamespace(
        components=[
            SimpleNamespace(feature_type="gene", feature_name="BAD", direction=None, weight=1.0),
            SimpleNamespace(feature_type="gene", feature_name="OK", direction=None, weight=1.0),
        ]
    )
    count, created = psc.create_signature_components_in_postgres(sig_model, sig, db=db)
    assert count == 1
    assert ("gene", "OK") in created
    assert db.commits == 1


def test_create_signature_components_outer_error_rolls_back(monkeypatch):
    db = FakeDB()
    sig_model = FakeSignatureModel(name="SigA")

    # Force db.query(SignatureComponentModel) to error by raising from query()
    class BoomDB(FakeDB):
        def query(self, model):
            raise RuntimeError("db boom")

    boom_db = BoomDB()
    monkeypatch.setattr(psc, "SignatureComponentModel", FakeComponentModel)

    sig = SimpleNamespace(components=[SimpleNamespace(feature_type="gene", feature_name="A", direction=None, weight=1.0)])
    count, created = psc.create_signature_components_in_postgres(sig_model, sig, db=boom_db)
    assert count == 0
    assert created == set()
    assert boom_db.rolled_back == 1


def test_create_signature_from_file_impl_calls_both(monkeypatch):
    db = FakeDB()
    sig = SimpleNamespace(name="SigA", description=None, modalities=[], components=[])

    created = FakeSignatureModel(name="SigA")
    monkeypatch.setattr(psc, "create_signature_in_postgres", lambda **k: created)
    monkeypatch.setattr(psc, "create_signature_components_in_postgres", lambda **k: (2, {("gene", "A")}))

    got = psc.create_signature_from_file_in_postgres(sig, db=db)
    assert got is created


def test_link_signature_to_dataset_success(monkeypatch):
    db = FakeDB()
    sig_id = uuid4()
    ds_id = uuid4()

    # Signature exists and dataset exists
    signature = FakeSignatureModel(name="SigA")
    signature.id = sig_id

    class FakeDataset:
        def __init__(self):
            self.signatures = []

    dataset = FakeDataset()

    # db.query(DatasetModel).first -> dataset; db.query(SignatureModel).first -> signature
    class DatasetModel:
        id = _Col("id")

    monkeypatch.setattr(psc, "SignatureModel", FakeSignatureModel)
    db.set_query(DatasetModel, first=dataset)
    db.set_query(FakeSignatureModel, first=signature)

    # Monkeypatch inside impl by injecting DatasetModel symbol via sys.modules replacement
    monkeypatch.setitem(sys.modules, "amprenta_rag.database.models", SimpleNamespace(Dataset=DatasetModel))

    ok = psc.link_signature_to_postgres_source(sig_id, source_type="dataset", source_id=ds_id, db=db)
    assert ok is True
    assert dataset.signatures == [signature]
    assert db.commits == 1


def test_link_signature_external_ids_fallback(monkeypatch):
    db = FakeDB()
    sig_id = uuid4()
    signature = FakeSignatureModel(name="SigA", external_ids={})
    signature.id = sig_id
    db.set_query(FakeSignatureModel, first=signature)

    monkeypatch.setattr(psc, "SignatureModel", FakeSignatureModel)

    ok = psc.link_signature_to_postgres_source(sig_id, source_type="dataset", source_notion_id="notion123", db=db)
    assert ok is True
    assert signature.external_ids and signature.external_ids["dataset_source"] == "notion123"
    assert db.commits == 1


def test_link_signature_exception_rollbacks(monkeypatch):
    class BoomDB(FakeDB):
        def query(self, model):
            raise RuntimeError("boom")

    db = BoomDB()
    monkeypatch.setattr(psc, "SignatureModel", FakeSignatureModel)
    ok = psc.link_signature_to_postgres_source(uuid4(), source_type="dataset", source_id=uuid4(), db=db)
    assert ok is False
    assert db.rolled_back == 1


