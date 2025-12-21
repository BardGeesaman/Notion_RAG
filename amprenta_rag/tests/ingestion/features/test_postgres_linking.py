from __future__ import annotations

from uuid import uuid4

from amprenta_rag.ingestion.features import postgres_linking as pl


class _FakeQuery:
    def __init__(self, obj=None, all_list=None):
        self.obj = obj
        self.all_list = all_list or []

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self.obj

    def all(self):
        return list(self.all_list)


class _FakeSession:
    def __init__(self, dataset=None, feature=None, literature=None, email=None):
        self.dataset = dataset
        self.feature = feature
        self.literature = literature
        self.email = email
        self.added = []
        self.committed = False

    def query(self, model):
        if model is pl.DatasetModel:
            return _FakeQuery(self.dataset)
        if model is pl.FeatureModel:
            return _FakeQuery(self.feature, [self.feature] if self.feature else [])
        if model is pl.Literature:
            return _FakeQuery(self.literature)
        if model is pl.Email:
            return _FakeQuery(self.email)
        return _FakeQuery()

    def add(self, obj):
        self.added.append(obj)

    def commit(self):
        self.committed = True

    def refresh(self, obj):
        return obj

    def rollback(self):
        self.committed = False


class _FakeFeature:
    def __init__(self, name, feature_type):
        self.name = name
        self.feature_type = feature_type
        self.id = uuid4()


class _FakeDataset:
    def __init__(self):
        self.features = []
        self.id = uuid4()


class _FakeRecord:
    def __init__(self):
        self.semantic_metadata = {}
        self.id = uuid4()


def test_link_features_to_dataset_appends(monkeypatch):
    dataset = _FakeDataset()
    db = _FakeSession(dataset=dataset, feature=_FakeFeature("f1", "metabolite"))
    monkeypatch.setattr(pl, "_find_or_create_feature_impl", lambda name, ftype, db: db.feature)

    pl.link_features_to_postgres_items(["F1"], dataset.id, "dataset", db=db)
    assert dataset.features  # feature appended


def test_link_features_updates_metadata(monkeypatch):
    record = _FakeRecord()
    db = _FakeSession(literature=record, feature=_FakeFeature("f1", "metabolite"))
    monkeypatch.setattr(pl, "_find_or_create_feature_impl", lambda name, ftype, db: db.feature)

    pl.link_features_to_postgres_items(["F1", "F2"], record.id, "literature", db=db)
    assert "feature_mentions" in record.semantic_metadata
    assert record.semantic_metadata["feature_mentions"]


def test_find_or_create_feature_creates_when_missing(monkeypatch):
    db = _FakeSession(feature=None)
    created = {}

    class FakeFeatureModel:
        def __init__(self, name, feature_type):
            created["name"] = name
            created["feature_type"] = feature_type
            self.name = name
            self.feature_type = feature_type
            self.id = uuid4()

    monkeypatch.setattr(pl, "FeatureModel", FakeFeatureModel)

    feature = pl.find_or_create_feature_in_postgres("New", "lipid", db=db)
    assert feature is not None
    assert created["name"] == "new"
    assert db.committed is True


def test_batch_link_features_to_dataset(monkeypatch):
    dataset = _FakeDataset()
    db = _FakeSession(dataset=dataset)

    class FakeFeatureModel:
        def __init__(self, name, feature_type):
            self.name = name
            self.feature_type = feature_type
            self.id = uuid4()

    monkeypatch.setattr(pl, "FeatureModel", FakeFeatureModel)

    results = pl.batch_link_features_to_dataset_in_postgres(
        features=[("feat1", "lipid"), ("feat2", "lipid")],
        dataset_id=dataset.id,
        db=db,
        max_workers=1,
    )
    assert all(results.values())
    assert len(dataset.features) == 2

